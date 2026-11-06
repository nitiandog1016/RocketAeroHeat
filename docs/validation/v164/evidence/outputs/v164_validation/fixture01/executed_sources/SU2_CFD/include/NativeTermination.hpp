#pragma once

// 模块用途。记录未处理的C++异常及Windows调用地址，随后保留原终止行为。
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <typeinfo>

namespace rah_termination {
inline std::terminate_handler previous_handler = nullptr;
inline volatile LONG recording = 0;

[[noreturn]] inline void record_and_terminate() noexcept {
  // 并发或递归终止时立即退出，避免异常处理再次进入求解器。
  if (InterlockedCompareExchange(&recording, 1, 0) != 0) std::_Exit(134);
  std::fprintf(stderr, "\nRAH_NATIVE_TERMINATION_BEGIN\nthread=%lu\n",
               static_cast<unsigned long>(GetCurrentThreadId()));
  const auto current = std::current_exception();
  if (current) {
    try {
      std::rethrow_exception(current);
    } catch (const std::exception& error) {
      std::fprintf(stderr, "exception_type=%s\nwhat=%s\n", typeid(error).name(), error.what());
    } catch (...) {
      std::fprintf(stderr, "exception_type=non-std\n");
    }
  } else {
    std::fprintf(stderr, "exception_type=none-current\n");
  }

  // 每帧同时给出模块基址和相对地址，供固定二进制符号表离线复核。
  void* frames[64] = {};
  const USHORT count = CaptureStackBackTrace(0, 64, frames, nullptr);
  for (USHORT index = 0; index < count; ++index) {
    HMODULE module = nullptr;
    char name[MAX_PATH] = {};
    const bool found = GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
        GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
        reinterpret_cast<LPCSTR>(frames[index]), &module) != 0;
    if (found) GetModuleFileNameA(module, name, MAX_PATH);
    const auto offset = reinterpret_cast<ULONG_PTR>(frames[index]) - reinterpret_cast<ULONG_PTR>(module);
    std::fprintf(stderr, "frame=%u address=%p module=%s base=%p offset=0x%llx\n",
                 static_cast<unsigned>(index), frames[index], name,
                 reinterpret_cast<void*>(module), static_cast<unsigned long long>(offset));
  }
  std::fprintf(stderr, "exit_action=delegate-original-terminate\nRAH_NATIVE_TERMINATION_END\n");
  std::fflush(stderr);
  if (previous_handler && previous_handler != record_and_terminate) previous_handler();
  std::abort();
}

inline void install() noexcept {
  previous_handler = std::set_terminate(record_and_terminate);
}
}  // namespace rah_termination

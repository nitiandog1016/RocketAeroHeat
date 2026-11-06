set pagination off
set confirm off
set print demangle on
break __cxa_call_terminate
break _ZN9__gnu_cxx24__trap_terminate_handlerEv
run
echo RAH_DEBUG_STOP\n
thread apply all bt 30
info registers rcx rdx r8 r9 rip
echo RAH_EXCEPTION_STORAGE\n
x/16gx $rcx-96
x/s *(char**)($rcx+40)
echo RAH_OBSERVATION_COMPLETE\n
quit

set pagination off
set confirm off
set print thread-events off
set environment _NO_DEBUG_HEAP 1
break __cxa_throw
commands
silent
printf "\nRAH_EXCEPTION\n"
x/s *(char**)($rdx+8)
bt 8
x/s *(char**)($rcx+16)
continue
end
run
thread apply all bt 12
quit

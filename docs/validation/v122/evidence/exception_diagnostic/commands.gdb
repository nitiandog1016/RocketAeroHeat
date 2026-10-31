set pagination off
set confirm off
set print thread-events off
break __cxa_throw
commands
silent
printf "\nRAH_EXCEPTION\n"
x/s *(char**)($rdx+8)
x/s *(char**)($rcx+8)
bt 8
continue
end
run
thread apply all bt 12
quit

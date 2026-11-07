import sys,time
for i in range(100):
 sys.stdout.write((str(i)+" observation row\n")*1000);sys.stdout.flush();time.sleep(.005)
print("SYNTHETIC-END",flush=True)

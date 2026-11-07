set pagination off
set confirm off
set print elements 8
set print max-depth 4
set debuginfod enabled off
break _D3lmr11fluidfvcell11FluidFVCell16signal_frequencyMFNibZd
run --max-cpus=1 --max-wall-clock=10
info args
python
import gdb,json,math
result={}
try:
    c=gdb.selected_frame().read_var('this').dereference()
    def scalar(x):
        y=float(x)
        return y if math.isfinite(y) else str(y)
    def vec(x): return [scalar(x[k]) for k in ('x','y','z')]
    def arr(x):
        names=[f.name for f in x.type.fields()]
        result.setdefault('sliceFields',names)
        return x['ptr'],int(x['length'])
    p,n=arr(c['lengths'])
    result['lengths']=[scalar(p[i]) for i in range(n)]
    result['Lmin']=scalar(c['L_min'])
    fs=c['fs'].dereference()
    result['velocity']=vec(fs['vel'])
    result['soundSpeed']=scalar(fs['gas']['a'])
    p,n=arr(c['iface'])
    result['normals']=[vec(p[i].dereference()['n']) for i in range(n)]
    result['selectedFaces']=[0,1,n-1]
    result['observed']=True
except Exception as e:
    result['error']=repr(e)
    print(gdb.execute('print *this',to_string=True))
open('time-observation.json','w').write(json.dumps(result,indent=2))
print('RAH-OBSERVATION',json.dumps(result))
end
finish
python
try:
    result['signal']=scalar(gdb.parse_and_eval('$xmm0.v2_double[0]'))
    open('time-observation.json','w').write(json.dumps(result,indent=2))
    print('RAH-SIGNAL',result['signal'])
except Exception as e: print('RAH-SIGNAL-ERROR',repr(e))
end
quit

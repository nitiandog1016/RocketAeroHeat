set pagination off
set confirm off
set debuginfod enabled off
break _D3lmr11fluidfvcell11FluidFVCell16signal_frequencyMFNibZd
run --max-cpus=1 --max-wall-clock=10
python
import gdb,json,struct,math
off=json.load(open('offsets.json'))
mem=gdb.selected_inferior()
def unpack(addr,fmt): return struct.unpack(fmt,mem.read_memory(addr,struct.calcsize(fmt)))
def value(x): return x if math.isfinite(x) else str(x)
def vector(addr): return list(unpack(addr,'<3d'))
c=int(gdb.parse_and_eval('$rdi'))
n,p=unpack(c+off['lengths'],'<2Q')
assert n==3
fs=unpack(c+off['fs'],'<Q')[0]
nf,pf=unpack(c+off['iface'],'<2Q')
result=dict(observed=True,abi=off,cellAddress=c,
    lengths=[value(x) for x in unpack(p,'<3d')],Lmin=value(unpack(c+off['L_min'],'<d')[0]),
    velocity=vector(fs+off['vel']),soundSpeed=unpack(fs+off['gas']+off['a'],'<d')[0],
    normals=[vector(unpack(pf+8*i,'<Q')[0]+off['normal']) for i in range(nf)],selectedFaces=[0,1,nf-1])
open('time-observation.json','w').write(json.dumps(result,indent=2))
print('RAH-OBSERVATION',json.dumps(result))
end
finish
python
result['signal']=value(float(gdb.parse_and_eval('$xmm0.v2_double[0]')))
open('time-observation.json','w').write(json.dumps(result,indent=2))
print('RAH-SIGNAL',result['signal'])
end
python
gdb.events.exited.connect(lambda e: open('inferior-exit.json','w').write(json.dumps(dict(exitCode=e.exit_code))))
end
continue
quit

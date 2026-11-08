// 模块：Time::run每步末及最终退出的只读状态审核；不回滚、不裁温、不改流动算子。
#include "functionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "rahPhaseDomain.H"
#include <cmath>
namespace Foam
{
namespace functionObjects
{
class rahPhaseStateGuard : public functionObject
{
    const fvMesh& mesh_;
    label lastIndex_,checks_;
    void check(const scalar p,const scalar T,const scalar rho,const scalar e,const vector& U,const label cell,const label patch) const
    {
        const scalar R=8.31451/.02896546;
        if (!rahPhaseDomain::accepted(p,T)||!std::isfinite(rho)||rho<=0||!std::isfinite(e)
            ||!std::isfinite(U.x())||!std::isfinite(U.y())||!std::isfinite(U.z())
            ||std::abs(rho/(p/(R*T))-1)>2e-10)
            FatalErrorInFunction<<"RAH_PHASE_FINAL_STATE_REJECT index="<<time_.timeIndex()
                <<" time="<<time_.value()<<" cellOrPatchFace="<<cell<<" patch="<<patch
                <<" p="<<p<<" T="<<T<<" rho="<<rho<<" e="<<e<<" U="<<U
                <<" pMax="<<rahPhaseDomain::ceiling(T)<<exit(FatalError);
    }
public:
    TypeName("rahPhaseStateGuard");
    rahPhaseStateGuard(const word& name,const Time& time,const dictionary& dict)
    : functionObject(name,time,dict),mesh_(time.lookupObject<fvMesh>(polyMesh::defaultRegion)),lastIndex_(-1),checks_(0)
    {
        if (dict.lookup<word>("model")!="rah-phase-aware-dry-air-of13-v1")
            FatalErrorInFunction<<"RAH_PHASE_GUARD_MODEL"<<exit(FatalError);
    }
    virtual wordList fields() const{return wordList{"p","T","rho","e","U"};}
    virtual bool execute()
    {
        const label index=time_.timeIndex();
        if (index==lastIndex_) return true;
        if (lastIndex_>=0&&index!=lastIndex_+1)
            FatalErrorInFunction<<"RAH_PHASE_GUARD_SKIPPED index="<<index<<" last="<<lastIndex_<<exit(FatalError);
        const auto& p=mesh_.lookupObject<volScalarField>("p");
        const auto& T=mesh_.lookupObject<volScalarField>("T");
        const auto& rho=mesh_.lookupObject<volScalarField>("rho");
        const auto& e=mesh_.lookupObject<volScalarField>("e");
        const auto& U=mesh_.lookupObject<volVectorField>("U");
        forAll(p,i) check(p[i],T[i],rho[i],e[i],U[i],i,-1);
        forAll(p.boundaryField(),patch)
            forAll(p.boundaryField()[patch],i)
                check(p.boundaryField()[patch][i],T.boundaryField()[patch][i],rho.boundaryField()[patch][i],
                    e.boundaryField()[patch][i],U.boundaryField()[patch][i],i,patch);
        lastIndex_=index;++checks_;
        if (index==0||time_.writeTime()||index%1000==0)
            Info<<"RAH_PHASE_GUARD_OK index="<<index<<" time="<<time_.value()<<" checks="<<checks_<<nl;
        return true;
    }
    virtual bool write(){return true;}
    virtual bool end()
    {
        execute();Info<<"RAH_PHASE_GUARD_END index="<<lastIndex_<<" checks="<<checks_<<nl;return true;
    }
};
defineTypeNameAndDebug(rahPhaseStateGuard,0);
addToRunTimeSelectionTable(functionObject,rahPhaseStateGuard,dictionary);
}
}

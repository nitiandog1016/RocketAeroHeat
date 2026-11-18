"""复核本轮端点修复，以及固定重放状态可复用的高精度误差数据。"""
import bisect
import csv
import math
import pathlib
from decimal import Decimal, localcontext

from runSU2RocketPhaseV222 import ROOT, read, save, sha
from auditCaloricInverseV260 import TABLE

BASE = ROOT / 'outputs/v260_validation'
OLD = BASE / 'case01'
CASE = BASE / 'case02'


def csv_rows(path):
    with path.open(newline='', encoding='utf8') as stream:
        return list(csv.DictReader(stream))


def as_decimal(text):
    return Decimal.from_float(float(text))


class TableIntervals:
    """在 80 位精度下积分每个表段，分别保存左、右端，显式识别表锚点舍入产生的缝隙和重叠。"""
    def __init__(self, kind, rows):
        gas_constant = as_decimal('287.0491267875601') if kind == 'E' else Decimal(0)
        self.knots = [row[0] for row in rows]
        self.lo, self.hi, self.c0, self.slope = [], [], [], []
        for (t, cp, h), (next_t, next_cp, _) in zip(rows, rows[1:]):
            width = next_t - t
            base = h - gas_constant*t + 300000
            capacity = cp - gas_constant
            assert capacity > 0 and next_cp - gas_constant > 0
            self.lo.append(base)
            self.hi.append(base + (capacity + next_cp - gas_constant)*width/2)
            self.c0.append(capacity)
            self.slope.append((next_cp - cp)/width)
        for values in (self.lo, self.hi):
            assert all(b > a for a, b in zip(values, values[1:]))

    def intervals(self, target):
        # 各段左闭右开，最高表温包含在末段内。
        first = bisect.bisect_left(self.hi, target)
        stop = bisect.bisect_right(self.lo, target)
        return [i for i in range(first, stop)
                if target < self.hi[i] or i == len(self.hi)-1]

    def roots(self, target, indices):
        temperatures = []
        for i in indices:
            lo, hi = self.knots[i:i+2]
            for _ in range(140):
                mid = (lo+hi)/2
                dt = mid - self.knots[i]
                value = self.lo[i] + self.c0[i]*dt + self.slope[i]*dt*dt/2
                if value < target:
                    lo = mid
                else:
                    hi = mid
            temperatures.append((lo+hi)/2)
        return temperatures


def audit():
    native = read(CASE/'native_completion.json')
    receipt = read(BASE/'process_edge-query/report.json')
    previous = read(BASE/'audit01/report.json')
    assert receipt['passed'] and receipt['inputsStable'] and not receipt['liveGroupAfterCleanup']
    assert native['normalReturn'] and native['samples'] == 195 and native['controlAcceptanceMismatches'] == 0
    assert native['entropyTemperatureBitMismatches'] == 0 and native['oldTracedMappingBitwiseOriginal22Fields']
    for name in ('states_input.csv', 'input_identity.json', 'states.csv', 'hs.csv'):
        assert sha(CASE/name) == previous['inputSHA256'][name] == sha(OLD/name), name

    inputs = csv_rows(CASE/'boundary_inputs.csv')
    tested = csv_rows(CASE/'boundary_results.csv')
    original = csv_rows(OLD/'controls.csv')
    updated = csv_rows(CASE/'controls.csv')
    assert len(original) == len(updated) == 60 and len(inputs) == len(tested) == 192
    assert len({r['id'] for r in inputs}) == len(inputs)
    comparison = []
    with localcontext() as ctx:
        ctx.prec = 80
        rows = [[as_decimal(x) for x in line.strip('{} ,').split(',')[:3]]
                for line in TABLE.read_text(encoding='utf8').splitlines()
                if line.startswith('{') and line[1].isdigit()]
        assert len(rows) == 3882
        tables = {kind: TableIntervals(kind, rows) for kind in ('E', 'H')}
        for source, observed in zip(inputs, tested):
            key = observed['id']
            target = float(observed['target'])
            assert key == source['id'] and observed['kind'] == source['kind']
            assert target.hex() == float(source['target']).hex()
            table = tables[observed['kind']]
            accepted = observed['newAccepted'] == '1'
            assert observed['newAccepted'] == observed['oldAccepted']
            assert observed['newReason'] == observed['oldReason']
            if key.startswith('original-'):
                i = int(key.split('-')[1])
                first, corrected = original[i], updated[i]
                assert observed['kind'] == first['kind'] and target.hex() == float(first['target']).hex()
                assert first['expectedValid'] == first['oldAccepted'] == observed['oldAccepted']
                assert corrected['newAccepted'] == corrected['oldAccepted'] == corrected['expectedValid']
                if first['label'] == 'node':
                    assert float(observed['newT']) == float(first['sourceT'])
            allowed = table.intervals(as_decimal(target)) if math.isfinite(target) else []
            entry = dict(id=key, accepted=accepted, mathematicalRootCount=len(allowed))
            if accepted:
                temperature = float(observed['newT'])
                assert 60 <= temperature <= 2000 and math.isfinite(temperature)
                if allowed:
                    roots = table.roots(as_decimal(target), allowed)
                    error = min(abs(as_decimal(temperature)-ref) for ref in roots)
                    assert error <= Decimal('2e-10'), (key, error)
                    entry['nearestRootErrorK'] = float(error)
                else:
                    nearest = min(table.knots, key=lambda t: abs(t-as_decimal(temperature)))
                    assert abs(nearest-as_decimal(temperature)) <= Decimal('2e-10'), key
                    entry['adjacentKnotK'] = float(nearest)
            else:
                assert math.isnan(float(observed['newT'])) and observed['newReason']
            comparison.append(entry)
        for kind in ('E','H'):
            for edge, valid_extreme in [('60', '+16'), ('2000', '-16')]:
                group = [r for r in tested if r['id'].startswith('band-'+kind+'-'+edge+'-')]
                assert len(group) == 33
                assert any(r['id'].endswith(valid_extreme) and r['newAccepted']=='1' for r in group)
                assert any(r['newAccepted']=='0' for r in group)
                assert [r['newAccepted'] for r in group] == sorted([r['newAccepted'] for r in group], reverse=edge=='2000')
        states = csv_rows(CASE/'states.csv')
        actual = [r for r in states if r['mode']=='new']
        assert all(len(tables['E'].intervals(as_decimal(r['receivedE'])))==1 for r in actual)
        hs = csv_rows(CASE/'hs.csv')
        for row in hs:
            assert abs(float(row['p'])/float(row['inputP'])-1) < 2e-10
            assert abs(float(row['s'])-float(row['targetS']))/max(1,abs(float(row['targetS']))) < 2e-10
            if row['mode'].endswith('H'):
                assert len(tables['H'].intervals(as_decimal(row['targetH']))) == 1

    recovered = tested[1]
    assert recovered['id']=='original-1' and original[1]['newAccepted']=='0' and recovered['newAccepted']=='1'
    assert previous['precisionPassed'] and previous['measurableGainPassed']
    assert previous['maximumErrorNewToOldRatio'] <= .5 and previous['pairedRMSErrorNewToOldRatio'] <= .5
    result = dict(passed=True, domainAcceptancePassed=True, caloricComponentAccepted=True,
                  sourceReplayRows=native['samples'], distinctReplayEnergyTargets=len({r['receivedE'] for r in actual}),
                  pairedDifferences=previous['pairedTemperatureDifferences'], frozenControls=60, extraEdgeNeighbors=132,
                  disagreements=0, fixedOriginallyRejectedRow=recovered, acceptedRowsWithoutUniqueMathematicalRoot=[
                      e['id'] for e in comparison if e['accepted'] and e['mathematicalRootCount']!=1],
                  maximumTemperatureErrorK=previous['maximumTemperatureErrorK'], pairedDeltaErrorRMSK=previous['pairedDeltaErrorRMSK'],
                  maximumErrorNewToOldRatio=previous['maximumErrorNewToOldRatio'],
                  pairedRMSErrorNewToOldRatio=previous['pairedRMSErrorNewToOldRatio'],
                  reusedByteIdenticalStatesAndHS=True, all195CaloricTargetsInUniqueInterval=True,
                  nativeSummary=native, priorAuditSHA256=sha(BASE/'audit01/report.json'),
                  CFDExecuted=False, newQualifiedAreaM2=0, fullSurfaceMilestoneComplete=False)
    out = BASE/'audit02'
    out.mkdir(exist_ok=False)
    save(out/'boundary_comparison.json',comparison)
    save(out/'report.json',result)
    print(result['fixedOriginallyRejectedRow'], result['disagreements'])


if __name__=='__main__':
    audit()

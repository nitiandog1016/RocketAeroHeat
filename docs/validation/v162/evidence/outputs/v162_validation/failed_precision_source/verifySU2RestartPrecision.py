"""独立验证弱热流附近双精度文本往返，保留旧保存态的严格失败身份。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from runEmbeddedSU2 import sha, atomic_write_json


def main():
    """模块一。用真实邻点温度和相邻ULP构造有区分力的序列化夹具。"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.relative_to(ROOT/'outputs')
    assert not output.exists()
    source = ROOT/'outputs/v161_validation/heat_discrepancy.json'
    state_path = ROOT/'outputs/v161_validation/control/restart_flow_000000.csv'
    row = json.loads(source.read_text())['top'][0]
    with state_path.open() as stream:
        reader = csv.reader(stream)
        fields = next(reader)
        for index, values in enumerate(reader):
            if index == row['neighbor']:
                original_token = values[fields.index('Temperature')].strip()
                break
    temperature = float(original_token)
    assert temperature == row['neighborTemperature']
    values = [temperature]
    for direction in (-np.inf, np.inf):
        value = temperature
        for _ in range(4):
            value = np.nextafter(value, direction)
            values.append(value)
    records = []
    for value in sorted(values):
        low, high = format(value, '.15e'), format(value, '.17g')
        low_read, high_read = float(low), float(high)
        factor = row['conductivity']/row['distance']
        records.append(dict(valueHex=float(value).hex(), scientificPrecision15=low,
            maxDigits10=high, roundTrip16=low_read == value, roundTrip17=high_read == value,
            lowPrecisionHeatChangeWm2=float((low_read-value)*factor),
            candidateNativeHeatDifferenceWm2=float((value-300)*factor-row['nativeHeat'])))
    assert all(r['roundTrip17'] for r in records)
    assert any(not r['roundTrip16'] for r in records)
    assert max(abs(r['lowPrecisionHeatChangeWm2']) for r in records) > 2e-10
    # 模块二。夹具说明机制，旧内存温度未保存，不能追溯认定旧热审核通过。
    report = dict(passed=True, originalTemperatureToken=original_token, sampleCount=len(records),
        roundTrip16Count=sum(r['roundTrip16'] for r in records), roundTrip17Count=len(records),
        records=records, originalMemoryTemperatureKnown=False, originalStrictHeatAuditPassed=False,
        physicsQualified=False, inputSHA256={str(p):sha(p) for p in (Path(__file__), source, state_path)})
    atomic_write_json(output, report)
    print(json.dumps({k:report[k] for k in ('passed', 'sampleCount', 'roundTrip16Count', 'roundTrip17Count')}))


if __name__ == '__main__':
    main()

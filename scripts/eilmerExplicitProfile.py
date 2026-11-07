"""显式局部伪时间方案：固定成熟算子、严格配置及最小接口准入。"""
PROFILE = 'rah-explicit-local-v1'
DEVELOPMENT_PROFILE = 'rah-explicit-local-development-v1'


def validate_explicit(value):
    """模块一：独立步数身份，不把显式更新称为Newton步或物理时间。"""
    if not isinstance(value, dict):
        raise ValueError('显式方案须为对象')
    development = value.get('profile') == DEVELOPMENT_PROFILE
    required = {'profile', 'maxSteps', 'checkpointStep'} | ({'candidate', 'candidateSHA256'} if development else set())
    if set(value) != required:
        raise ValueError('显式方案字段不完整或含未知参数')
    if value['profile'] not in (PROFILE, DEVELOPMENT_PROFILE):
        raise ValueError('未支持的显式方案')
    for key in ('maxSteps', 'checkpointStep'):
        if type(value[key]) is not int or value[key] <= 0:
            raise ValueError('显式步数须为正整数')
    if development:
        import re
        if value['maxSteps'] != 200 or value['checkpointStep'] != 100:
            raise ValueError('当前具名发展仅支持登记的200步/100检查点')
        if not isinstance(value['candidate'],str) or not value['candidate'].startswith('/') or not re.fullmatch('[0-9a-f]{64}',value['candidateSHA256']):
            raise ValueError('完整域缺少准确候选引用')
    elif not value['checkpointStep'] < value['maxSteps'] <= 8:
        raise ValueError('当前仅开放最小显式接口；完整域/长任务须另登记准入')
    return dict(value)


def settings(value):
    """模块二：完整列出改变的推进/输出参数，保留原有物理与空间离散。"""
    value = validate_explicit(value)
    return dict(solver_mode='time-marching', gasdynamic_update_scheme='tvd_rk3',
                with_local_time_stepping=True, local_time_stepping_limit_factor=1000,
                with_super_time_stepping=False, residual_smoothing=False,
                stringent_cfl=True, cfl_value=0.1, cfl_count=1,
                viscous_signal_factor=1.0, turbulent_signal_factor=1.0,
                fixed_time_step=False, dt_init=1e-10, dt_max=1e-3,
                max_time=1.0, max_step=value['maxSteps'],
                max_attempts_for_step=1, print_count=1, write_transient_residuals=True,
                write_loads=True, dt_plot_value=1.0, dt_loads=1.0,
                write_loads_at_step=value['checkpointStep'],
                write_flow_solution_at_step=value['checkpointStep'])


def make_explicit_job(nk_job, condition):
    """模块三：只复用已经审核的物理/网格Lua段，移除全部NK配置。"""
    import json
    prefix = nk_job.split('NewtonKrylovGlobalConfig', 1)[0]
    prefix = prefix.replace('config.solver_mode = "steady"', 'config.solver_mode = "time-marching"')
    prefix += '-- 显式TVD RK3与局部伪时间；此时间不代表真实轨迹时间。\n'
    for key, value in settings(condition['explicit']).items():
        prefix += 'config.%s = %s\n' % (key, json.dumps(value))
    return prefix


def audit_explicit_config(cfg, condition, control):
    """模块四：审核原生实际输出；Lua拼接正确本身不能证明配置生效。"""
    expected = settings(condition['explicit'])
    merged = dict(cfg, **control)
    # 固定版本将推进控制分存control，并将标量周期写成调度数组。
    for key, value in expected.items():
        if key in ('dt_plot_value', 'cfl_value'):
            continue
        if merged.get(key) != value:
            raise ValueError('显式原生配置不符: ' + key)
    for prefix, value in [('cfl', 0.1), ('dt_plot', 1.0)]:
        if any(cfg.get(prefix + suffix) != target for suffix, target in
               [('_schedule_length', 1), ('_schedule_times', [0.0]), ('_schedule_values', [value])]):
            raise ValueError('显式CFL/保存调度不符')
    return True

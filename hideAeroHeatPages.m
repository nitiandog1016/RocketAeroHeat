function hideAeroHeatPages(host)
%HIDEAEROHEATPAGES 切换页面时隐藏旧控件，窗口关闭时统一释放。
% MATLAB 浏览器 UI 的布局事件可能延迟到达；立即删除旧控件会导致
% 已删除对象的布局回调错误。隐藏后保留其生命周期至窗口关闭。
% 同时禁用旧页面的交互，防止旧工况控件在后台更改结果。
pages = host.Children;
for i=1:numel(pages)
    pages(i).Visible = 'off';
    controls = findall(pages(i));
    for k=1:numel(controls)
        if isprop(controls(k),'Enable'), controls(k).Enable='off'; end
    end
end
drawnow;
end

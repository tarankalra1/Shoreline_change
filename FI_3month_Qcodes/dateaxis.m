function dateaxis(str)
datetick('x',str,'keepticks')
h = zoom;
set(h,'ActionPostCallback',{@zoom_date,str});
h_pan = pan;
set(h_pan,'ActionPostCallback',{@zoom_date,str});
end

function zoom_date(obj,event_obj,str)

figlist=get(gcf,'Children');
for i=1:length(figlist)
    if (strcmpi(get(figlist(i),'type'),'axes'))&(~strcmpi(get(figlist(i),'tag'),'legend'))
        break; % find the last axis in current figure
    end
end
set(figlist(i),'xTickMode','Auto','xticklabelMode','Auto')
datetick(figlist(i),'x',str,'keepticks')
end


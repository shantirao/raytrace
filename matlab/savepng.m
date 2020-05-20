function filename = savepng( name, xscale, yscale, size)
%savepng Save a .PNG image of the current figure
set(gcf,'PaperPositionMode', 'auto')
set(gcf,'PaperType','usletter')
oldPosition = get(gcf,'Position');
if nargin < 4
    size = [1200 675];
end
if nargin > 1
    size(1) = size(1) * xscale;
    if nargin > 2
        size(2) = size(2) * yscale;
    else
        size(2) = size(2) * xscale;
    end
end

set(gcf, 'Position', [128 128 size]); %[128 128 1280 1024]);
filename = [name '.png'];

%fill second monitor
if exist('\Program Files\ImageMagick-6.8.1-Q16\convert.exe','file')
print(gcf,'-r150','-dpng','print.png');
!"\Program Files\ImageMagick-6.8.1-Q16\convert.exe" -trim print.png print.png
copyfile('print.png',filename);
elseif exist('e:\Program Files\ImageMagick-6.9.1-Q16\convert.exe','file')
print(gcf,'-r150','-dpng','print.png');
!"e:\Program Files\ImageMagick-6.9.1-Q16\convert.exe" -trim print.png print.png
copyfile('print.png',filename);
else
print(gcf,'-r150','-dpng',filename);
end

set(gcf,'Position', oldPosition);
end


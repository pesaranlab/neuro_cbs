function imhandle = tvimage(image, varargin)
%  TVIMAGE displays images in a better way than imagesc
%
%  TVIMAGE(IMAGE) displays the IMAGE in the proper orientation
%  and allows you to edit the labels and ranges with a series
%  of optional arguments:
%
%  'XLabel' - String giving the text to label the X-axis
%
%  'YLabel' - String giving the text to label the Y-Axis
% 
%  'XNTicks' - Integer giving the number of X tickmarks
%
%  'YNTicks' - Integer giving the number of Y tickmarks
% 
%  'XRange' - Two element vector giving the minmax for the X-axis
% 
%  'YRange' - Two element vector giving the minmax for the Y-axis
%
%  'CLim' - Two element vector giving the minmax for the color scale
%
%   For example, to define the color scale:
%       tvimage(IMAGE,'CLim',[0,0.4])
%
%   or to define both the color scale and the XRange
%       tvimage(IMAGE,'CLim',[0,0.4],'XRange',[0,100])
%

%default values
XLabel='';
YLabel='';
XNTicks=5;
YNTicks=5;

image=squeeze(image);

sIm=size(image);
XRange=[0,sIm(1)];
YRange=[0,sIm(2)];
CLim=[];

assign(varargin{:});

xticklabel=sprintf('%2.1f',XRange(1));
for x = 1:XNTicks
    xticklabel=[xticklabel,sprintf('|%2.1f',...
            x*(XRange(2)-XRange(1))./XNTicks+XRange(1))];
end

yticklabel=sprintf('%2.1f',YRange(2));
for y=YNTicks-1:-1:0
    yticklabel=[yticklabel,sprintf('|%2.1f',...
            y*(YRange(2)-YRange(1))./YNTicks+YRange(1))];
end

xtickpos=1;
if XNTicks > 0
    for x=1:XNTicks
        xtickpos=[xtickpos,x*sIm(1)./XNTicks];
    end
end

if YNTicks > 0
    ytickpos=1;
    for y=1:YNTicks
        ytickpos=[ytickpos,y*sIm(2)./YNTicks];
    end
end

NX = size(image,1);
NY = size(image,2);
X = linspace(XRange(1),XRange(2),NX);
Y = linspace(YRange(1),YRange(2),NY);

if size(image,3) > 1
    imhandle = imagesc(X,Y,permute(image,[2 1 3]));
    axis xy
elseif length(CLim) > 0
    imhandle = imagesc(X,Y,image',CLim);
    axis xy
else
    imhandle = imagesc(X,Y,image');
    axis xy
end
% set(gca,'XTick',xtickpos);
% set(gca,'YTick',ytickpos);     
set(get(gca,'Ylabel'),'String',YLabel);
set(get(gca,'Xlabel'),'String',XLabel);
% set(gca,'YTickLabel',yticklabel)
% set(gca,'XTickLabel',xticklabel)



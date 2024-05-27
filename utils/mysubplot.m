function h=mysubplot(N,M,n)

figure('DefaultAxesXTickLabel','','DefaultAxesYTickLabel','');
set(gcf,'WindowState','maximized');
defPos = get(gcf,'DefaultAxesPosition');
XOffset = 0.01;
YOffset = 0.01;
Left = defPos(1);
Bottom = defPos(2);
Width = (defPos(3)-(M-1)*XOffset)/M;
Height = (defPos(4)-(N-1)*YOffset)/N;
axLeft = Left : Width+XOffset: Left + (M-1)*(Width+XOffset);
axBottom = Bottom : Height+YOffset : Bottom + (N-1)*(Height+YOffset);
[axLeft,axBottom] = meshgrid(axLeft,axBottom);
h = zeros(N*M,1);
for i = 1:N*M
   h(i) = axes('Position',[axLeft(i) axBottom(i) Width Height]);
end
if nargin>2
    axes(h(n));
end

for i=10139:10139
    I=imread(['',num2str(i),'.tif']);
    C=imagesc(I);
    [A, ret]=imcrop(C);
    E=numel(A);
    F(i)=sum(A(:))./E;

end
for i=10000:10139
    I=imread(['',num2str(i),'.tif']);
    B=imcrop(I,ret);
    F(i)=sum(B(:))./E;

end
G=F(10000:10139);
G=G-min(G);

% H=G./max(G);
% MH=max(H)/2;
% 
% 
% for i=1:139
%     if (H(i+1)-MH >= 0 && H(i)-MH <= 0)
%     break;  
%     end  
% end
% locate=i 
% G=smooth(G);
% dG=diff(G);
% G=G./max(G);
filename = 'int4.xlsx';
xlswrite(filename,(G'));
plot(0:5/60:139*5/60,G)
ylabel('\bfIntensity[a.u.]','FontSize',14);
xlabel('\bfTime[min]','FontSize',14);
set(gca,'Fontsize',14,'FontWeight','bold','linewidth',1.5);

% filename = 'normalize.xlsx';
% xlswrite(filename,(H'));
% plot(0:5/60:139*5/60,smooth(H))
% ylabel('\bfIntensity[a.u.]','FontSize',14);
% xlabel('\bfTime[min]','FontSize',14);
% set(gca,'Fontsize',14,'FontWeight','bold','linewidth',1.5);
% figure 
% filename = 'speed4.xlsx';
% xlswrite(filename,(dG'));
% plot(smooth(dG))
% ylabel('\bfThe rate of intensity[a.u.]','FontSize',14);
% xlabel('\bfFrame','FontSize',14);
% % set(gca,'Fontsize',9,'FontWeight','bold','linewidth',1.5);

%%
%  J�Ǹ���excel�е�6��sheet�����ݣ���BBB����dx�У����Ҳ���û�����в������ĵ�
k=linspace(0.01,0.99,99);
p=linspace(0.01,0.99,99);
S=[S5(1:99,5)'];
for i=1:1:98
    A=S5((1+i*99):(99+i*99),5)';
    S=[S;A];
end
figure(1);
%subplot(2,2,3),
surf(k,p,S);
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('sw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);


F=[S5(1:99,4)'];
for i=1:1:98
    A=S5((1+i*99):(99+i*99),4)';
    F=[F;A];
end
figure(2);
%subplot(2,2,3),
surf(k,p,F);
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('fw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);
%%
%�õ���Ȼ������Ĳ����������ڸ���excel�е�7��sheet����BBB������Ч�ú����У��в������ĵ�
%��B����������
k=linspace(0.01,0.99,99);
p=linspace(0.01,0.99,99);
S=[S6(1:99,5)'];
for i=1:1:98
    A=S6((1+i*99):(99+i*99),5)';
    S=[S;A];
end
figure(3);
%subplot(2,2,3),
surf(k,p,S);
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('sw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);


F=[S6(1:99,4)'];
for i=1:1:98
    A=S6((1+i*99):(99+i*99),4)';
    F=[F;A];
end
figure(4);
%subplot(2,2,3),
surf(k,p,F);
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('fw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);
%%
%��8��sheet�������˲��������н�����������ģ�BBB������Ч�ú�����
k=linspace(0.01,0.99,99);
p=linspace(0.01,0.99,99);
S=[S1(1:99,5)'];
for i=1:1:98
    A=S1((1+i*99):(99+i*99),5)';
    S=[S;A];
end
figure(5);
%subplot(2,2,3),
surf(k,p,S');
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('sw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);


F=[S1(1:99,4)'];
for i=1:1:98
    A=S1((1+i*99):(99+i*99),4)';
    F=[F;A];
end
figure(6);
%subplot(2,2,3),
surf(k,p,F');
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('sw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);
p=1;
for j=2:98
    for i=2:98
        if ( S(i,j)>S(i+1,j) && S(i,j)>S(i-1,j) )&&( S(i,j)>S(i,j+1) && S(i,j)>S(i,j-1))
            jjj(p,1)=i;
            jjj(p,2)=j;
            p=p+1;
        end
    end
end
% %%
% a=1;
% b=1;
% c=1;
% d=1;
% if (a==1 && b==1) &&  (c==1 && d==1)
%     disp('lll')
% else
%     disp('kkk')
% end
%%
k=linspace(0.01,0.99,99);
p=linspace(0.01,0.99,99);
S=[S7(1:99,5)'];
for i=1:1:98
    A=S7((1+i*99):(99+i*99),5)';
    S=[S;A];
end
figure(7);
%subplot(2,2,3),
surf(k,p,S');
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('sw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);


F=[S7(1:99,4)'];
for i=1:1:98
    A=S7((1+i*99):(99+i*99),4)';
    F=[F;A];
end
figure(8);
%subplot(2,2,3),
surf(k,p,F');
%grid off % ȥ���������� 
shading interp % ȥ��ͼ���ϵ����񣬼�ʹ֮�⻬
view([45 35]);
xlabel('p','Fontsize',14);
ylabel('k','Fontsize',14);
zlabel('sw','Fontsize',14);
xlim([0.01 0.99]);
ylim([0.01 0.99]);

p=1;
for j=2:98
    for i=2:98
        if ( S(i,j)>S(i+1,j) && S(i,j)>S(i-1,j) )&&( S(i,j)>S(i,j+1) && S(i,j)>S(i,j-1))
            ppp(p,1)=i/100;
            ppp(p,2)=j/100;
            p=p+1;
        end
    end
end
clear all; 
close all;
%clc;

tic;
count=1;
for k=0.41:0.01:0.99
    for p=0.01:0.01:0.99
        BBB=0        ;
        disp(k);
        disp(p);
        improper=0;
        for cc=1:3
            gamma=0.3; %-3会使n和V变得很平   %-2
            epsilon=-0.1;%好像不能太大  %-0.05  %epsilon的绝对值稍微大一些，可以使n上升的幅度变大，V也会有一些改变

            rho=0.5;%rho变小会使n的变化幅度变小
            C=0.1;%C越小，断点越大，n的幅度也大   %2
            % k=0.99%k越小，断点越大   %0.9
            % p=0.99                  %0.4  

            a=0.8;
            eta=0.6;
            xmean=a/eta;
            sigma=0.8;
            sd=sigma^2/2/eta;
            xmin=xmean-2*sd;
            xmax=xmean+2*sd;

            I=1000;
            hmin=0.0001;  %-0.02
            hmax=100;
            h=linspace(hmin,hmax,I)';   %列向量
            dh=(hmax-hmin)/(I-1);

            J=7;
            x=linspace(xmin,xmax,J);    %行向量
            dx=(xmax-xmin)/(J-1);
            dx2=dx^2;

            hh=h*ones(1,J);
            xx=ones(I,1)*x;

            mu=a-eta.*x;%+BBB;
            s2=0.5*sigma.^2;

            maxit=100;
            crit=10^(-6);
            Delta=1000;

            Vhf = zeros(I,J);
            Vhb = zeros(I,J);
            Vxf = zeros(I,J);
            Vxb = zeros(I,J);
            Vxx = zeros(I,J);
            n = zeros(I,J);

            %CONSTRUCT MATRIX Bswitch SUMMARIZING EVOLUTION OF x
            chi =  - min(mu,0)/dx + s2/dx2;              
            upsilon =  min(mu,0)/dx - max(mu,0)/dx - s2*2/dx2;
            zeta = max(mu,0)/dx + s2/dx2;

            %This will be the upperdiagonal of the B_switch    B的上对角线
            updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
            for j=1:J
                updiag=[updiag;repmat(zeta(j),I,1)];%repmat把zeta(j)复制I行1列
            end

            %This will be the center diagonal of the B_switch    B的中间对角线
            centdiag=repmat(chi(1)+upsilon(1),I,1);
            for j=2:J-1
                centdiag=[centdiag;repmat(upsilon(j),I,1)];
            end
            centdiag=[centdiag;repmat(upsilon(J)+zeta(J),I,1)];

            %This will be the lower diagonal of the B_switch
            lowdiag=repmat(chi(2),I,1);
            for j=3:J
                lowdiag=[lowdiag;repmat(chi(j),I,1)];
            end

            %Add up the upper, center, and lower diagonal into a sparse matrix
            Bswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

            v0=(xx.*hh.^(1-gamma)/(1-gamma))/rho;%-C/(1+1/epsilon);
            v=v0;

            maxit=30;
            for m=1:maxit
                V=v;
                Vhf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/dh;
                Vhf(I,:)=(V(I,:)-V(I-1,:))/dh;    %boundary
               %Vhf(I,:) = C/(1-k).*(p/(1-k).*hmax).^(1/epsilon);

                Vhb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/dh;
                Vhb(1,:) =(V(2,:)-V(1,:))/dh;   % = 0; %boundary   %我可不可以从这个作差里求出nf呢？
                %Vhb(1,:) = C/(1-k).*(1/(1-k)+(p/(1-k)).*hmin).^(1/epsilon);

                I_concave = Vhb > Vhf ;%indicator whether value function is concave (problems arise if this is not the case)


                 %consumption and savings with forward difference
                nf=((1-k).*Vhf/C).^epsilon;
                hf=(1-k).*nf-p.*hh;%这里加上了一个+1

                % backward difference
                nb=((1-k)*Vhb/C).^epsilon;
                hb=(1-k).*nb-p.*hh;%这里加上了一个+1

                %稳态时n和价值函数的导的值
                n0=p/(1-k).*hh;%这里加上了一个   -1/(1-k)
                Vh0=C/(1-k).*n0.^(1/epsilon);

                % dV_upwind makes a choice of forward or backward differences based on
                % the sign of the drift   
                If=hf>0;
                Ib=hb<0;
                I0=(1-If-Ib);


                Vh_Upwind=Vhf.*If + Vhb.*Ib + Vh0.*I0;%important to include third term 

                n=((1-k)/C.*Vh_Upwind).^epsilon;%这个有没有改进的机会。
                u=BBB+xx.*hh.^(1-gamma)/(1-gamma)-C/(1+1/epsilon).*n.^(1+1/epsilon); 

                X = - min(hb,0)/dh;
                %X(1,:)=0;
                Y = - max(hf,0)/dh + min(hb,0)/dh;
                Z = max(hf,0)/dh;
                %Z(I,:)=0;

                updiag=0; %This is needed because of the peculiarity of spdiags.
                for j=1:J
            %         disp('updiag=');
                    updiag=[updiag;Z(1:I-1,j);0];    
                end
                %updiag;

                centdiag=reshape(Y,I*J,1) ;   %按照列来reshape

                lowdiag=X(2:I,1);
                for j=2:J  %从2到10闭区间
                    lowdiag=[lowdiag;0;X(2:I,j)];
                end

                %稀疏矩阵
                AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);

                A = AA + Bswitch;

            %     disp(max(abs(sum(A,2))));
                if max(abs(sum(A,2)))>10^(-9)%对行求和得到的列向量  %这个是不是判断收敛
                   disp('Improper Transition Matrix');
                   improper=1;
                   break
                end

                B = (1/Delta + rho)*speye(I*J) - A;%B是I*J行，I*J列的稀疏方矩阵

                u_stacked = reshape(u,I*J,1);%u本来是I行J列的矩阵，reshape为I*J行1列的列向量
                V_stacked = reshape(V,I*J,1);

                b = u_stacked + V_stacked/Delta;%b为I*J行1列的列向量

                V_stacked = B\b; %左除：%SOLVE SYSTEM OF EQUATIONS

                %V_stacked2=inv(B)*b;

                V = reshape(V_stacked,I,J);

                Vchange = V - v;
                v = V;

                dist(m) = max(max(abs(Vchange)));    %第二次。第一个max返回包含了每列最大值的行向量，第二个max返回最大值。dist是一个向量而不是一个函数
            %     disp(dist(m));
%                 disp(m);
                if dist(m)<crit
                    disp('Value Function Converged, Iteration = ')
                    disp('')
                    disp(m)

                    break
                end
            end
            toc;    

            % hd = (1-k).*n-p.*hh;      %hd是对h(t)求导
            % set(gca,'FontSize',14)
            % figure(1)
            % %subplot(2,2,1),
            % plot(h,hd,h,zeros(1,I),'--');
            % xlim([hmin hmax]);
            % xlabel('h')
            % ylabel('hd(h,x)')
            % % % print -depsc derivativehabit.eps
            % % 
            % set(gca,'FontSize',14)
            % figure(2)
            % %subplot(2,2,2),
            % plot(h,n,h,zeros(1,I),'--')
            % xlim([hmin  hmax])
            % xlabel('h')
            % ylabel('n(h,x)')
            % % % print -depsc stationarytime2D.eps
            % % 
            % % X=real(X);
            % % Y=real(Y);
            % % Z=real(Z);
            % % x=real(x);
            % % 
            % % icut=990;
            % % hcut=h(1:icut);
            % % ncut=n(1:icut,:);
            % % 
            % % hcut=real(hcut);
            % % ncut=real(ncut);
            % % 
            % % figure(3);
            % % %subplot(2,2,3),
            % % surf(hcut,x,ncut');
            % % %grid off % 去掉坐标网格 
            % % shading interp % 去掉图像上的网格，即使之光滑
            % % view([45 35]);
            % % xlabel('Habit, h','Fontsize',14);
            % % ylabel('Environment, x','Fontsize',14);
            % % zlabel('Payment times, n(h.x)','Fontsize',14);
            % % xlim([hmin hmax]);
            % % ylim([xmin xmax])
            % % % print -depsc stationarytime3D.eps
            % % 
            % figure(4);
            % set(gca,'FontSize',14)
            % %subplot(2,2,4),
            % plot(h,V,h,zeros(1,I),'--');
            % xlabel('h')
            % ylabel('V(h,x)')
            % xlim([hmin hmax]);
            % % % print -depsc valuefunction.eps
            % figure(6);
            % set(gca,'FontSize',14)
            % %subplot(2,2,4),
            % plot(h,u,h,zeros(1,I),'--');
            % xlabel('h')
            % ylabel('u(h,x)')
            % xlim([hmin hmax]);

            AT=A';
            b=zeros(I*J,1);

            i_fix=5;
            b(i_fix)=0.1;
            row=[zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
            AT(i_fix,:)=row;

            gg=AT\b;
            g_sum=gg'*ones(I*J,1)*dh*dx;
            %%%
            gg=gg./g_sum;

            g=reshape(gg,I,J);

            % figure(5);
            % shading interp;
            % icut=100;
            % hcut=h(1:icut);
            % gcut=g(1:icut,:);
            % surf(hcut,x,gcut');
            % view([45 35])
            % xlabel('Habit,h','FontSize',14)
            % ylabel('Environment,x','FontSize',14)
            % zlabel('Density g(h,x)','FontSize',14)
            % xlim([hmin max(hcut)])
            % ylim([xmin xmax])
            % alpha(.7)

            % print -depsc densities_diffusion.eps
            % %%
            % disp('AT')
            % disp(AT(1:5,1:5))
            % disp(AT(6995:7000,6995:7000))

            T=k.*n+p.*hh;
            % T(find(isnan(T)==1)) = 0;
            TotalFirmWelfare=ones(1,I)*(T.*g)*ones(J,1);
            TotalSocialWelfare=ones(1,I)*(V.*g)*ones(J,1);
            if BBB~=TotalFirmWelfare
                BBB=TotalFirmWelfare;
            else
                break;
            end
        end
        Record(count,1)=k;
        Record(count,2)=p;
        Record(count,3)=BBB;
        Record(count,4)=TotalFirmWelfare;
        Record(count,5)=TotalSocialWelfare;
        Record(count,6)=improper;
        disp(count);
        count=count+1;
        
    end
end
 %%%
% disp('a') ;
% disp(2)
 % count=1;
% for k=0:1:10
%     for p=20:1:30
%        J(count,1)=k;
%        J(count,2)=p;
%        count=count+1; 
%        if p==25
%            J(count,2)='improper';
%        end
%     end
% end
% if 2~=1
%     disp('e')
% end
% %%
% k=J(:,1);
% p=J(:,2)
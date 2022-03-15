function [ainf,avac1,avac2,Fpop,acase,ahosp,adth]=world_odes(pop0,dems,Ieff1,Ieff2,Seff1,Seff2,Deff1,Deff2,VD1,VD2,numwaned,case_risk,control,hscal,dscal,sev_risk,wane_symp,wane_trans,wane_sev,dth_risk,varscale,age_mat,tau,MaxTime)

%reduction in susceptibility due to vaccination
Ivacc=1-(Ieff1*(VD1-VD2)+Ieff2*(VD2-numwaned'));
Ivacc16=[Ivacc(1:15)',sum(Ivacc(16:21)'.*(dems(16:21)'/sum(dems(16:21))))]';

%reduction in probability of symptoms due to vaccination
CR=case_risk.*((1-VD1')+Seff1.*(VD1'-VD2')*(1-Ieff1)+Seff2.*(VD2'-wane_symp*numwaned)*(1-Ieff2))./((1-VD1')+(VD1'-VD2')*(1-Ieff1)+(VD2'-wane_symp*numwaned)*(1-Ieff2));
CR16=[CR(1:15),sum(CR(16:21).*dems(16:21)'/sum(dems(16:21)))];
%reduction in transmissibility due to vaccination
transad=(1-VD1')+(VD1'-wane_trans*numwaned)*0.55;
transad=[transad(1:15),sum(transad(15:end).*dems(15:21)'/sum(dems(15:21)))];

%transmission matrix
agemat=diag(varscale*Ivacc16')*age_mat*diag(CR16+tau*(1-CR16))*diag(transad)*min(max(control,0.1),10);
agemat(:,end)=agemat(:,end)*sum(dems(16:21))/dems(16);
agemat=reshape(agemat,1,[]);

gamma=1;
lag=1;
m=5;
N=[dems(1:15);sum(dems(16:21))];
MaxTime=max(MaxTime,1);

%solve ODEs
options = odeset('RelTol', 1e-5,'NonNegative',1);
[~, pop]=ode45(@Diffs,[0:1:MaxTime],pop0,options,[agemat,N',gamma,lag,m]);


%collect results
ainf=zeros(size(pop,1),21);
ainf(:,1:15)=repmat(dems(1:15)',size(pop,1),1)-pop(:,1:15);
binf=ones(size(pop,1),1)*N(16);
binf=repmat((binf-pop(:,16))./binf,1,21-16+1);
ainf(:,16:21)=binf.*repmat(dems(16:21)',size(binf,1),1);
ainf(ainf<0)=0;

avac1=repmat(VD1',MaxTime,1);
avac2=repmat(VD2',MaxTime,1);

Sev=sev_risk.*((1-VD1')+Deff1.*(VD1'-VD2')*(1-Ieff1)*Seff1+Deff2.*(VD2'-wane_sev*numwaned)*(1-Ieff2)*Seff2)./((1-VD1')+(VD1'-VD2')*(1-Ieff1)*Seff1+(VD2'-wane_sev*numwaned)*(1-Ieff2)*Seff1);

acase=(ainf(2:end,:)-ainf(1:(end-1),:)).*repmat(CR,size(ainf,1)-1,1);
ahosp=acase.*repmat(Sev*hscal,size(ainf,1)-1,1);
adth=acase.*repmat(Sev.*dth_risk*hscal*dscal,size(ainf,1)-1,1);

Fpop=pop(end,:);

%ODES
    function dpop=Diffs(t, pop, parameter)
        
        L=16;
        a=1:L;
        
        agemat=reshape(parameter(1:(L^2)),L,L);
        N=parameter((L^2+1):(L^2+L))';
        gamma=parameter(L^2+L+1);
        lag=parameter(L^2+L+2);
        m=parameter(L^2+L+3);
        
        S=zeros(L,1);I=S;
        E=zeros(L,m);
        S(a)=pop(a);
        for i=1:m
            E(a,i)=pop(i*L+a);
        end
        I(a)=pop((m+1)*L+a);
        
        dpop=0*pop;
        dpop(a)=-agemat*(I./N).*S;
        dpop(a+L)=agemat*(I./N).*S-E(:,1)*lag;
        for i=2:m
            dpop(a+i*L)=E(:,i-1)*lag-E(:,i)*lag;
        end
        dpop(a+(m+1)*L)=E(:,m)*lag-gamma*I;
        
        
    end


end

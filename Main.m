

%%
clear all

load globalDATA
load globalPARAMS
load wanePARAMS
load Contacts



%% Choose strategies to simulate 

thresholds=[100,65,40,15]; %sharing thresholds (100=full sharing)
adapted=1; %adapted behaviour on=1, off=0


%% Other parameters

COUNTRIES=1:length(Country); %which countries to simulate
whovacc=[0,0,3/5,1,ones(1,17)];%5 year age groups to vaccinate
uptake=[0.8*ones(1,12),0.9*ones(1,9)];%uptake 90% above 60, 80% below
past=365; %days after 01/01/20 to start simulation
T=datenum(2022,1,1)-datenum(2020,1,1)-past; %end point for simulation
int=2; %simulation step size
trange=floor(T/int); %number of simulation steps
its=1; %number of simulation samples
clag=2;% control smoothing for adapted behaviour

%% arrays to store results
INF=zeros(length(Country),T+past,length(thresholds));
HOS=zeros(length(Country),T+past+5,length(thresholds));
DTH=zeros(length(Country),T+past+5,length(thresholds));
CAS=zeros(length(Country),T+past,length(thresholds));
VAC1=zeros(T,21,length(Country),length(thresholds));
VAC2=zeros(T,21,length(Country),length(thresholds));

 %%  Run simulation
for vs=1:length(thresholds)

     
    %initialise values
    V_orig=V_origs;
    vaccleft=vacclefts;
    boostleft=0*DEM1;
    VD1=VD1s;
    VD2=VD2s;
    cDone=Inf.*ones(1,length(Country));
    NUMWANED=NUMWANEDs;
    varscale=VARIANT(2,find(VARIANT(1,:)>sum(INFs(:,end),'all')/sum(DEM1,'all'),1,'first'));
    VARSCALE=ones(1,trange)*varscale;
    CONTROL=CONTROLs;
    CONTROL2=CONTROL; %reference for adapted behaviour
    
    I=round(INF_strt'*1e-3); 
    E=repmat(I,1,10);
    S=[DEM1(1:15,:);sum(DEM1(16:21,:))]'-INF_strt(:,:)'*(1+0.01/10);
    POP=[S,E,I]';
    
    
    %set scale of fitted control values for each country for adapted behaviour
    ctscale=zeros(15,length(Country));
    for i=COUNTRIES
        if min(CONTROL(1:183,i))==max(CONTROL(1:183,i))
            ctscale(:,i)=1;
        else
            ctscale(:,i)=min(CONTROL(1:183,i)):(max(CONTROL(1:183,i))-min(CONTROL(1:183,i)))/14:max(CONTROL(1:183,i));
        end
    end
    
    %redistributed vaccine
    V_share=zeros(size(V_orig,1),1);
    
    inf=zeros(T,21,length(Country));
    hos=zeros(T,21,length(Country));
    cas=zeros(T,21,length(Country));
    dth=zeros(T,21,length(Country));
    vac1=zeros(T,21,length(Country));
    vac2=zeros(T,21,length(Country));
    propvacc=1- vaccleft./(uptake'.*2.*DEM1);
    propvacc(isnan(propvacc))=0;
    
    %sharing on/off 1/0
    ctswitch=zeros(1,length(Country));
    %control scale index for behaviour adapting
    ctval=ones(1,length(Country));
    
    
    for t_int=1:trange %simulation time interval
        t_day=t_int*int+past; %days from start of 2020
        %                 t/trange
        fprintf(1,'%0.1f%% , threshold=%d  \n', round(100*t_int/trange,1),thresholds(vs));
        
        
        for i= COUNTRIES
            if ~isnan(missing(i))
                
                
                %update vaccine still needed
                vaccleft(:,i)=max(whovacc'.*DEM1(:,i).*uptake'.*2*0.99-(VD1(:,i)+VD2(:,i)).*DEM1(:,i),0);
                %proportion of each agegroup vaccinated
                propvacc_new=min(1- vaccleft(:,i)./(uptake'.*2.*DEM1(:,i)),1);   propvacc_new(isnan(propvacc_new))=0;
                %newly vaccinated proportion in interval
                diff_pvacc=propvacc_new-propvacc(:,i);
                propvacc(:,i)=propvacc_new;
                
                %add waning of newly vaccinated
                NUMWANED(t_day+1:end,:,i)= NUMWANED(t_day+1:end,:,i)+cumsum(diff_pvacc'.*[wane_profile(1:min(size(NUMWANED,1)-t_day,wane_time)),zeros(1,size(NUMWANED,1)-t_day-min(size(NUMWANED,1)-t_day,wane_time))]');
                
                %find waned and eligable for boosters
                NUMWANED_boost=NUMWANED;
                NUMWANED_boost(NUMWANED_boost./repmat(uptake,size(NUMWANED,1),1,size(NUMWANED,3))<wane_am*min(180/wane_time,1))=0; %only boost after 6 months
                
                %country vaccine need
                vneed=(1-propvacc)+squeeze(NUMWANED_boost(t_day,:,:))./uptake';
                
                
                vneed(vneed<0)=0;
                if sum(vneed(:,i))~=0
                    %vaccines given = nonshared vaccine + shared vaccine allocation
                    V_sim=sum(V_orig(past+t_int:past+t_int+int-1,i))+sum(V_share(past+t_int:past+t_int+int-1))*sum(vneed(:,i).*DEM1(:,i))/(sum(vneed(:,:).*DEM1(:,:),'all'));
                else
                    %if vaccination completed, all future is shared
                    V_sim=0;
                    
                    V_share=V_share+V_orig(:,i);
                    V_orig(:,i)=0;
                    
                end
                
                %identify individuals to be given booster vaccine
                if sum(vaccleft(:,i))<V_sim && cDone(i)==Inf
                    cDone(i)=t_day;
                    boostleft(:,i)=DEM1(:,i).*whovacc'.*uptake';
                end
                
                
                
                
                deliv_boost=0;
                if V_sim~=0 % if country has vaccine, find who gets vaccinated
                    num2vacc=cumsum(flipud(vaccleft(:,i)));
                    v_indx=find(num2vacc>V_sim,1,'first'); % index of youngest age group that gets vaccinated in this time interval
                    if isempty(v_indx)
                        deliv_boost=V_sim-num2vacc(end);
                        V_sim=1;v_indx=21;
                    elseif v_indx==1
                        V_sim=[zeros(21-v_indx,1);V_sim/(num2vacc(v_indx));ones(v_indx-1,1)];
                    else
                        V_sim=[zeros(21-v_indx,1);(V_sim-num2vacc(v_indx-1))/(num2vacc(v_indx)-num2vacc(v_indx-1));ones(v_indx-1,1)];
                    end
                    
                    %use remaining vaccine allocation to boost
                    num2vacc=cumsum(flipud(boostleft(:,i)));
                    v_indx_b=find(num2vacc>deliv_boost,1,'first');
                    if isempty(v_indx_b)
                        deliv_boost=ones(21,1);
                    elseif v_indx_b==1
                        deliv_boost=[zeros(21-v_indx_b,1);deliv_boost/(num2vacc(v_indx_b));ones(v_indx_b-1,1)];
                    else
                        deliv_boost=[zeros(21-v_indx_b,1);(deliv_boost-num2vacc(v_indx_b-1))/(num2vacc(v_indx_b)-num2vacc(v_indx_b-1));ones(v_indx_b-1,1)];
                    end
                    boostleft(:,i)=boostleft(:,i).*(1-deliv_boost);
                    
                    
                    if v_indx>21-(thresholds(vs)/5+1) && ctswitch(i)==0
                        V_share=V_share+V_orig(:,i);
                        V_orig(:,i)=0;
                        ctswitch(i)=1;
                        [~,ctval1]=min(abs(ctscale(:,i)-CONTROL(t_int,i)));
                        ctval(i)=ctval1(1);
                    end
                    
                end
                
                
                
                %update who has vaccine in each country
                VD1(:,i)=VD1(:,i)+max(whovacc'.*uptake'-VD1(:,i),0).*V_sim;
                VD2(:,i)=VD2(:,i)+max(whovacc'.*uptake'-VD2(:,i),0).*V_sim;
                
                %boosted no longer waned
                NUMWANED(t_day:end,:,i)=(1-deliv_boost').*NUMWANED(t_day:end,:,i)+deliv_boost'.*[wane_profile(1:min(size(NUMWANED,1)-t_day+1,wane_time)),zeros(1,size(NUMWANED,1)-t_day+1-min(size(NUMWANED,1)-t_day+1,wane_time))]';
                
                %run ODES
                [ainf,avac1,avac2,pop0,acase,ahosp,adth]=world_odes(POP(:,i),DEM1(:,i),Ieff1,Ieff2,Seff1,Seff2,Deff1,Deff2,VD1(:,i),VD2(:,i),NUMWANED(t_day,:,i),case_risk,mean(CONTROL(t_int:t_int+clag-1,i)),hscal(i),dscal(i),sev_risk,wane_symp,wane_trans,wane_sev,dth_risk,VARSCALE(t_int),diag(sigma16(i,:))*K(:,:,i),tau,int);
                
                pop0=max(pop0,0);
               
                if adapted==1
                    if ctswitch(i)==1 %adapted behaviour if sharing
                        if sum(pop0((1:16*5)+16))>1.02*sum(POP((1:16*5)+16,i)) %increasing infection
                            ctval(i)=max(ctval(i)-1,1); %relax control
                        elseif sum(pop0((1:16*5)+16))*1.02<sum(POP((1:16*5)+16,i)) %decreasing infection
                            [~,ctval1]=min(abs(ctscale(:,i)-CONTROL2(t_int+clag,i))); %only allow control above current
%                               ctval(i)=min(ctval(i)+1,ctval1);
                                   ctval(i)=min(ctval(i)+1,size(ctscale,1));
                        end
                        CONTROL(t_int+clag,i)=ctscale(ctval(i),i); %update future control value
                    end
                end
                
                %collect values
                 POP(:,i)=pop0;
                inf([1:int]+(t_int-1)*int,:,i)=ainf(2:end,:);
                vac1([1:int]+(t_int-1)*int,:,i)=avac1;
                vac2([1:int]+(t_int-1)*int,:,i)=avac2;
                
                cas([1:int]+(t_int-1)*int,:,i)=acase;
                hos([1:int]+(t_int-1)*int,:,i)=ahosp;
                dth([1:int]+(t_int-1)*int,:,i)=adth;
            end
            
            
        end
        %update variant transmission scaling
        varscale=VARIANT(2,find(VARIANT(1,:)>sum(inf(t_int*int,:,:),'all')/sum(DEM1,'all'),1,'first'));
        VARSCALE(t_int:end)=varscale;
        
    end
    
%collect results
    
    INF(:,:,vs)=[INFs,squeeze(sum(inf(:,:,:),2))'];
    HOS(:,:,vs)=[zeros(179,9),[HOSs,(HOSs(:,past)+cumsum(squeeze(sum(hos(1:end-4,:,:),2)))')]];
    DTH(:,:,vs)=[zeros(179,20),[DTHs,(DTHs(:,past)+cumsum(squeeze(sum(dth(1:end-15,:,:),2)))')]];
    CAS(:,:,vs)=[CASs,(CASs(:,past)+cumsum(squeeze(sum(cas,2)))')];
    VAC1(:,:,:,vs)=vac1;
    VAC2(:,:,:,vs)=vac2;
    
end

%% Save results in table

Results=tblerun(INF,DTH,VAC2,Country,DEM1,inc_brk,thresholds);

writetable(Results,'Results.csv')

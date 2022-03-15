function [Results]=tblerun(INF,DTH,VAC2,country,DEM1,inc_brk,thresholds)
%% Function to store results in table

rdate=datenum(2022,1,1)-datenum(2020,1,1);

gCountry=find(sum(INF(:,:,1),2)~=0);
 
  Country=string(repelem(country(gCountry),length(thresholds))');
  
  inclabs={'Low income','Lower middle income','Upper middle income','High income'};
   Income_Bracket=strings(length(gCountry),1);
  for i=1:length(gCountry)
  Income_Bracket(i)=inclabs{inc_brk(gCountry(i))};
  end
  Income_Bracket=repelem(Income_Bracket,length(thresholds));

  Sharing_strategy=strings(1,length(thresholds));
  for i=1:length(thresholds)
   Sharing_strategy(i)=sprintf('%d+ threshold' , thresholds(i));
  end
   Sharing_strategy=repmat(Sharing_strategy,1,length(gCountry))';
   
   Proportion_Vaccinated=reshape(squeeze(sum(squeeze(VAC2(rdate-365,:,gCountry,:)).*DEM1(:,gCountry),1))'./sum(DEM1(:,gCountry),1),length(gCountry)*length(thresholds),1);
   Proportion_Infected=reshape(squeeze(INF(gCountry,rdate,:))'./sum(DEM1(:,gCountry),1),length(gCountry)*length(thresholds),1);
   Deaths_per_100000=reshape(squeeze(DTH(gCountry,rdate,:))'./sum(DEM1(:,gCountry),1)*1e5,length(gCountry)*length(thresholds),1);
   
 
    Results=table(Country,Income_Bracket,Sharing_strategy,Proportion_Vaccinated, Proportion_Infected,Deaths_per_100000);
    
    

function [spgm,ys]=centerSpectrogramF(spgm,ntl,nLenses,ys,halfShift,thresh)
% Center spectrogram
% if halfshift==1, will shift by an additional half lens
spgmOffset=spgm-min(min(spgm));
fProj=sum(spgmOffset,2); fProj=fProj-min(fProj);
rise=find(fProj>0.6*max(fProj),1); fall=find(fProj>0.6*max(fProj),1,'last');

% fProj(fProj<0.6*max(fProj))=0;
fProj(1:rise)=0; fProj(fall:end)=0;
centerF=round((1:ntl)*fProj/sum(fProj));
totalShift=0;
shift=round(ntl/2-centerF-halfShift*round(ntl/2));
totalShift=totalShift+shift;

while abs(shift)>thresh
ys=circshift(ys,shift);
spgm=reshape(ys,ntl,nLenses);
fProj=sum(spgm,2);fProj(fProj<0.6*max(fProj))=0;
centerF=round((1:ntl)*fProj/sum(fProj));
totalShift=0;
shift=round(ntl/2-centerF);
totalShift=totalShift+shift;

end

% % % 
% % % % Center spectrogram
% % % % if halfshift==1, will shift by an additional half lens
% % % % if the amount to shift that is found is smaller then tresh (i.e., an
% % % % index) then the program will terminate
% % % spgmOffset=spgm-min(min(spgm));
% % % fProj=sum(spgmOffset,2); fProj=fProj-min(fProj);
% % % rise=find(fProj>0.2*max(fProj),1); fall=find(fProj>0.2*max(fProj),1,'last');
% % % 
% % % % fProj(fProj<0.6*max(fProj))=0;
% % % fProj(1:rise)=0; fProj(fall:end)=0;
% % % centerF=round((1:ntl)*fProj/sum(fProj));
% % % totalShift=0;
% % % shift=round(ntl/2-centerF-halfShift*round(ntl/2));
% % % totalShift=totalShift+shift;
% % % % 
% % % % while abs(shift)>thresh
% % % ys=circshift(ys,shift);
% % % spgm=reshape(ys,ntl,nLenses);
% % % % fProj=sum(spgmOffset,2);fProj(fProj<0.6*max(fProj))=0;
% % % % centerF=round((1:ntl)*fProj/sum(fProj));
% % % % shift=round(ntl/2-centerF-halfShift*round(ntl/2));
% % % % totalShift=totalShift+shift;
% % % % 
% % % % end

end

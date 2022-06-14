
function [spgm,ys]=centerSpectrogramF(spgm,ntl,nLenses,ys,halfShift,thresh)
% Center spectrogram
% if halfshift==1, will shift by an additional half lens
fProj=sum(spgm,2);fProj(fProj<0.6*max(fProj))=0;
centerF=round((1:ntl)*fProj/sum(fProj));
totalShift=0;
shift=ntl/2-centerF-halfShift*round(ntl/2);
totalShift=totalShift+shift;

while abs(shift)>thresh
ys=circshift(ys,shift);
spgm=reshape(ys,ntl,nLenses);
fProj=sum(spgm,2);fProj(fProj<0.6*max(fProj))=0;
centerF=round((1:ntl)*fProj/sum(fProj));
totalShift=0;
shift=ntl/2-centerF;
totalShift=totalShift+shift;

end


end

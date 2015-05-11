function new_img=labeling(A,B,npix)
new_img=A;
% first: aortic valve
M = A.*B + A;
M=uint8(M);
bound=false(size(M));
%second: PVs

nvox=3;
delta=nvox-1;

bound(1:(1+delta),:,:)=true;
bound((end-delta):end,:,:)=true;
bound(:,1:(1+delta),:)=true;
bound(:,(end-delta):end,:)=true;
bound(:,:,1:(1+delta))=true;
bound(:,:,(end-delta):end)=true;


M=M+ uint8(2*(M & bound));


new_img=uint8(M);


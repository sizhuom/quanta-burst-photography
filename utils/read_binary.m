function frames=read_binary(file_path,file_name,bin_suffix,rows,cols)

file_size=8*1024*1024*8; %8MB
no_frames=file_size/rows/cols;

fileID = fopen(file_path+ "/" + string(floor(bin_suffix/2^11)) + "/" + string(floor(bin_suffix/2^6)) + "/" + file_name + string(bin_suffix) + ".bin");
% tic
A = fread(fileID);
B = uint8(dec2bin(A))-48;
fclose(fileID);
% 
% Ba=reshape(B,[8388608 4 2]);
% Bb=fliplr(Ba);
% Bc=Bb(:,:,[2 1]);
% Bd=reshape(Bc,[8388608 8]);
% Be=reshape(Bd,[8388608*8 1]);
Bp = reshape(permute(B,[2 1]),[file_size 1]);
C=reshape(Bp,[no_frames,rows,cols]);
% D=permute(C,[3 2 1]);
E=reshape(C, [4, 2,no_frames*rows*cols/8]);
F=flip(E);
G=F(:,[2 1],:);
H=reshape(G,[cols rows no_frames]);
HH=permute(H,[3 2 1]);
% HHH=reshape(HH,[67108864/rows/cols rows cols]);
frames=HH;
% toc
% figure; image(squeeze(sum(uint16(HH),1)));

% B=reshape(A,[512, 256, 512]);
% C=permute(B,[3 2 1]);
% frames=C;

% f = waitbar(0);
% for i=1:min(no_frames,frames)
% %     i
%     waitbar(i/no_frames,f);
%     for j=1:rows
%         for k=1:(cols/8)
%             tmp=B(k+(j-1)*(cols/8)+(i-1)*(cols/8)*rows,:);            
%             frames(i,j,((k-1)*8+1):((k-1)*8+4))=fliplr(tmp(5:8));            
%             frames(i,j,((k-1)*8+5):((k-1)*8+8))=fliplr(tmp(1:4));            
%         end
%     end
% end
% 
% close(f)

end
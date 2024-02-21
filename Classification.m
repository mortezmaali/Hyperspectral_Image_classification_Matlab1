datafile = 'C:\Users\Morteza\OneDrive\Desktop\PhD\New_Data\8cal_Seurat_AFTER';
hdrfile = 'C:\Users\Morteza\OneDrive\Desktop\PhD\New_Data\8cal_Seurat_AFTER.hdr';

hcube = hypercube(datafile,hdrfile);

img = im2double(hcube.DataCube(:,:,1:end-1));
[m n l] = size(img);

sig = zeros(7,l);
a = img(46,164,:);
sig(1,:) = reshape(a,[1 l]) ;


a = img(389,677,:);
sig(2,:) = reshape(a,[1 l]) ;

a = img(316,174,:);
sig(3,:) = reshape(a,[1 l]) ;

a = img(68,309,:);
sig(4,:) = reshape(a,[1 l]) ;

a = img(46,92,:);
sig(5,:) = reshape(a,[1 l]) ;

a = img(66,1051,:);
sig(6,:) = reshape(a,[1 l]) ;

a = img(565,72,:);
sig(7,:) = reshape(a,[1 l]) ;

hcuben = hypercube(img,hcube.Wavelength(1:end-1));

abundanceMap = estimateAbundanceLS(hcuben,sig','Method','fcls');

classNames = {'class 1';'class 2';'class 3';'class 4';'class 5';'class 6';'class 7'};

[~,matchIdx] = max(abundanceMap,[],3);
figure
imagesc(matchIdx)
colorbar('Ticks',1:7,'TickLabels',classNames)
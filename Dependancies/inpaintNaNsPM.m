function PM2 = inpaintNaNsPM(PM)

a = 1*exp(1i*PM);
PM2 = zeros(size(PM));
for i=1:size(PM,3)
    PM2(:,:,i) = angle(inpaint_nans(a(:,:,i),3));
end



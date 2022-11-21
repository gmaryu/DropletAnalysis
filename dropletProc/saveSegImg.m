function saveSegImg(i, imnc, imcy)
nameFormat1='img_000000%03d_nuc_000.jpg';
fname1= sprintf(nameFormat1,i-1);
nameFormat1='img_000000%03d_cyto_000.jpg';
fname2= sprintf(nameFormat1,i-1);
imwrite(imnc, fname1);
imwrite(imcy, fname2);
end
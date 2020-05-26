[x,y] = meshgrid((-10:10),(-10:10));

c = [0 0 0 0 1];

[z dxdy] = zernikeXY(c, [x(:) y(:)], 10);

subplot(1,4,1);
imagesc(reshape(z,21,21),[-1 1]); axis image; colorbar;

subplot(1,4,2);
imagesc(reshape(dxdy(:,1),21,21),[-1 1]); axis image; colorbar;


subplot(1,4,3);
imagesc(reshape(dxdy(:,2),21,21),[-1 1]); axis image; colorbar;

subplot(1,4,4);

imagesc(reshape(dxdy(:,1).^2 + dxdy(:,2).^2,21,21),[-1 1]/10);axis image; colorbar;
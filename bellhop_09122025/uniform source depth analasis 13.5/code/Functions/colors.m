function res = colors(filename)
img = double(imread(filename));
img = reshape(img, [length(img) 3]);
ind = round(linspace(1, length(img), 256));
res = img(ind, :) / 255;
end


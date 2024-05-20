%% noiseonf function credit:
% Copyright (c) 1996-2009 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%
% Full noiseonf available in MATLAB link

noiseonf([1550 1550], 1.11);
ans = ans - min(ans(:));
ans = ans/max(ans(:));ans;
mean = mean(ans, "all");
ans = im2bw(ans, mean);
imwrite(ans, "image.jpg");



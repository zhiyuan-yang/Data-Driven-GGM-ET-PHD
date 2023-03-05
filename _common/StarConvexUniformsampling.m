function y = StarConvexUniformsampling(X,Rot,sizeObject,R,N)

%measurement cov: R
%sampling number: N

a = sizeObject(1);
b = sizeObject(2);
c = sizeObject(3);
d = sizeObject(4);

y = zeros(size(X,1),N);
for s = 1:N
    s_x = -a + 2 * a * rand();
    s_y = -c + 2 * c * rand();
    notInside = s_y > b && s_x < -d  ||...
               s_y > b && s_x > d   ||...
               s_y < -b && s_x < -d ||...
               s_y < -b && s_x > d ; 
                
    while(notInside)
        s_x = -a + 2 * a * rand();
        s_y = -c + 2 * c * rand();
        notInside = s_y > b && s_x < -d  ||...
                   s_y > b && s_x > d   ||...
                   s_y < -b && s_x < -d ||...
                   s_y < -b && s_x > d ; 
    end
    y(:,s) = Rot*[s_x;s_y] + X + (randn(1,2) * chol(R))';
end

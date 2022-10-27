function plotter(xi,pos,q,colo)
x = pos(3:3:end);
y = pos(1:3:end);

xhi = xi(2:3:end);

l = 1;
if nargin>2
    l = q;
end
col = 'b';
if nargin>3
    col = colo;
end
for node = 1:length(x)
    end1 = [0;l/2];
    end2 = [0;-l/2];
    M = [cos(xhi(node)) -sin(xhi(node));
        sin(xhi(node)) cos(xhi(node))];
    end1_n = [x(node);y(node)] + M*end1;
    end2_n = [x(node);y(node)] + M*end2;
    plot([end1_n(1),end2_n(1)],[end1_n(2),end2_n(2)],col);
    hold on;
end
end
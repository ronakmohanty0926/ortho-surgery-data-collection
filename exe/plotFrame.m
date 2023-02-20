function plotFrame(p1,p2,p3)


line([p1(1,1),0],[p1(1,2),0],[p1(1,3),0],'Color','red');
hold on;
line([p2(1,1),0],[p2(1,2),0],[p2(1,3),0], 'Color','green');
line([p3(1,1),0],[p3(1,2),0],[p3(1,3),0], 'Color','blue');
daspect([1,1,1]);

function plotFrames(data)

figure;
view(3);

line([1,0],[0,0],[0,0],'Color','red','LineWidth',2);
hold on;
line([0,0],[1,0],[0,0], 'Color','green', 'LineWidth',2);
line([0,0],[0,0],[1,0], 'Color','blue', 'LineWidth',2);
daspect([1,1,1]);

xlabel('X');
ylabel('Y');
zlabel('Z');

% for i = 3:length(data(:,1))
for i = 1:1
    
p1 = data(i,1:3);
p2 = data(i+1,1:3);
p3 = data(i+2,1:3);

plotFrame(p1, p2, p3)
 
end
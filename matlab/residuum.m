function residuum
[meshPath, simulPath, outputPath] = getPath;
residuum = load([outputPath,'residuum.txt']);

figure('color','w')
box on
plot(residuum(:,1),(residuum(:,3)-residuum(1,3))/residuum(1,3)*100,'linewidth',2)
xlabel('iter')
ylabel('total mass change [%]')
grid on

figure('color','w')
box on
plot(residuum(:,1),residuum(:,2),'linewidth',2)
xlabel('iter')
ylabel('u_{max} [m/s]')
grid on
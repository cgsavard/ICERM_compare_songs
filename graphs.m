

figure(2)
imagesc(PIs{6,1})

figure(3)
birth = songs{6,1}(:,1);
death = songs{6,1}(:,2);
hold on
title('Length vs. Midpoint')
xlabel('Midpoint');
ylabel('Length');
scatter((death + birth)/2,(death - birth));

figure(4)
hold on
title('Length vs. Start')
xlabel('Midpoint');
ylabel('Length');
scatter(birth, death-birth);


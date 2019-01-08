
a = 5
figure(5)
imagesc(PIs{a,1})

figure(6)
birth = songs{a,1}(:,1);
death = songs{a,1}(:,2);
hold on
title('Length vs. Midpoint')
xlabel('Midpoint');
ylabel('Length');
scatter(((death + birth)/2)/(max(death)-min(birth)), (death - birth));

figure(7)
hold on
title('Length vs. Start')
xlabel('Start');
ylabel('Length');
scatter(birth/(max(birth)-min(birth)), death-birth);


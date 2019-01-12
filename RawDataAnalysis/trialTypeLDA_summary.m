% slow data
mean_ = [0.7189    0.8829    0.9458;
         0.7164    0.9405    0.9707;
         0.6529    0.9188    0.9393;
         0.6867    0.9405    0.9979;
         0.6312    0.8781    0.9791;
        ];


std_  = [0.0861    0.0520    0.0426;
         0.1570    0.0329    0.0111;
         0.1681    0.0239    0.0309;
         0.1119    0.0418    0.0035;
         0.0775    0.0722    0.0091;
        ];

ngroups = 3;
nbars   = 5;

figure
hold on
bar(1:ngroups, mean_', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_(i,:), std_(i, :), 'k', 'linestyle', 'none');
end
set(gca, 'TickDir', 'out')
xlim([0.5, 3.5])
ylim([0.5 1.05])
set(gca, 'XTick', 1:3)
set(gca, 'YTick', 0.5:0.1:1)
setPrint(4*3, 3, 'decode_bar_slow', 'pdf')
    

    
% fast data    
mean_ = [0.6248    0.6547    0.7119;
         0.5480    0.6095    0.7857;
         0.5249    0.5846    0.7761;
        ];


std_  = [0.0476    0.0513    0.0512;
         0.0394    0.0513    0.1012;
         0.0297    0.0498    0.0593;
        ];

ngroups = 3;
nbars   = 3;

figure
hold on
bar(1:ngroups, mean_', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_(i,:), std_(i, :), 'k', 'linestyle', 'none');
end
set(gca, 'TickDir', 'out')
xlim([0.5, 3.5])
ylim([0.5 1.05])
set(gca, 'XTick', 1:3)
set(gca, 'YTick', 0.5:0.1:1)
setPrint(4*3, 3, 'decode_bar_fast', 'pdf')
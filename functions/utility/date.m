function str = date()

c = clock;
d = sprintf('%02.0f',c(3)+100);
m = sprintf('%02.0f',c(2));
y = sprintf('%02.0f',c(1));
str = [y(3:4) '_' m '_' d(2:3)];

end


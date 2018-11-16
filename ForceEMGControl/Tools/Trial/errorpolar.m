function[] = errorpolar(angles,meann,stdd,color)
w = 0.05;
for i = 1:length(angles)
    polarplot([angles(i) angles(i)],[meann(i)-stdd(i),meann(i)+stdd(i)],color)
    polarplot([angles(i)-w angles(i)+w],[meann(i)-stdd(i),meann(i)-stdd(i)],color)
    polarplot([angles(i)-w angles(i)+w],[meann(i)+stdd(i),meann(i)+stdd(i)],color)
end
end
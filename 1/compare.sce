clear;
clc;

function ds = gravity(t, s)
    mu = 398600.4418;
    ds = [s(4:6); -mu*s(1:3)/norm(s(1:3))^3];
endfunction

mu = 398600.4418;
a = 7000;
n = sqrt(mu/a^3);
period = 2*%pi/n;
t = linspace(0, 2*period, 721);
th = n*t;
xset = [0.1 10 100];
yset = [0 0.1 1];
zset = [0 0.1];
root = get_absolute_file_path("compare.sce");
out = fullfile(root, "results");
if ~isdir(out) then
    mkdir(out);
end
f = mopen(fullfile(out, "comparison.csv"), "wt");
mfprintf(f, "case,variant,x0_km,y0_km,z0_km,vx0_km_s,vy0_km_s,max_error_km,final_error_km,max_separation_km,ratio\n");
// Reference spacecraft initial state in an inertial frame whose axes coincide with RTN at t=0.
ref0 = [a;0;0;0;sqrt(mu/a);0];
ref = ode(ref0, 0, t, 1d-10, 1d-10, gravity);
caseid = 0;
for ix = 1:3
    for iy = 1:3
        for iz = 1:2
            x0 = xset(ix);
            y0 = yset(iy);
            z0 = zset(iz);
            for mode = 1:5
                caseid = caseid+1;
                if mode == 1 then
                    vx0 = n*y0/2;
                    vy0 = -2*n*x0;
                    variant = "Согласованная";
                else
                    vxvals = [0 0 0.1 0.1];
                    vyvals = [0 0.1 0 0.1];
                    vx0 = vxvals(mode-1);
                    vy0 = vyvals(mode-1);
                    variant = "Несогласованная";
                end
                // v_inertial = v_reference + v_relative + omega cross rho.
                dep0 = ref0+[x0;y0;z0;vx0-n*y0;vy0+n*x0;0];
                dep = ode(dep0, 0, t, 1d-10, 1d-10, gravity);
                rel = zeros(3, length(t));
                for j = 1:length(t)
                    R = ref(1:3, j)/norm(ref(1:3, j));
                    H = cross(ref(1:3, j), ref(4:6, j));
                    N = H/norm(H);
                    T = cross(N, R);
                    rel(:, j) = [R T N]'*(dep(1:3, j)-ref(1:3, j));
                end
                x = (4-3*cos(th))*x0+sin(th)/n*vx0+2*(1-cos(th))/n*vy0;
                y = y0+6*(sin(th)-th)*x0-2*(1-cos(th))/n*vx0+(4*sin(th)-3*th)/n*vy0;
                z = z0*cos(th);
                hcw = [x;y;z];
                err = sqrt(sum((hcw-rel).^2, "r"));
                sep = sqrt(sum(rel.^2, "r"));
                mfprintf(f, "%d,%d,%.3f,%.3f,%.3f,%.10f,%.10f,%.9f,%.9f,%.9f,%.9f\n", caseid, mode, x0, y0, z0, vx0, vy0, max(err), err($), max(sep), max(sep)/a);
                figure(1);
                clf();
                gcf().figure_size = [1100 850];
                gcf().background = color("white");
                subplot(2, 2, 1);
                plot(rel(1, :), rel(2, :), "b-");
                plot(x, y, "r--");
                xgrid();
                titleText = [msprintf("Случай %d: %s скорость", caseid, variant);
                    msprintf("x0=%.1f км, y0=%.1f км, z0=%.1f км", x0, y0, z0);
                    msprintf("vx0=%.3f м/с, vy0=%.3f м/с, vz0=0 м/с", vx0*1000, vy0*1000)];
                xtitle(titleText, "x, км", "y, км");
                legend(["Кеплер" "ХКУ"]);
                subplot(2, 2, 2);
                plot(t/period, err);
                xgrid();
                xtitle("Расхождение положений", "Время, витки", "Ошибка, км");
                subplot(2, 2, 3);
                plot(t/period, rel(3, :), "b-");
                plot(t/period, z, "r--");
                xgrid();
                xtitle("Поперечное движение", "Время, витки", "z, км");
                subplot(2, 2, 4);
                plot(t/period, sep);
                xgrid();
                xtitle("Относительное расстояние", "Время, витки", "Расстояние, км");
                xs2pdf(gcf(), fullfile(out, msprintf("case_%02d.pdf", caseid)));
            end
        end
    end
end
mclose(f);
mprintf("Готово: %d вариантов, 2 витка, результаты в %s\n", caseid, out);

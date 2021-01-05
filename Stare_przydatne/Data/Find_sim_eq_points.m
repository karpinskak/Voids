for k=1:numel(promien)
                    [par_eq]=wylicz_param(Const,par_set,promien(k),part(1).par.delta,part(1).par.teta,part(1).par.A,1,1);
                    [points, stability] = eq_points(part(1).par.A,par_eq.Sv,par_eq.St);
                    punkty=[punkty;points];
                    stabilnosc=[stabilnosc,stability];
end

punkty_x=skala*punkty(:,1)*part(1).par.delta;
punkty_y=skala*punkty(:,2)*part(1).par.delta;
scatter(punkty_x,punkty_y,10);

                    
                    
                    
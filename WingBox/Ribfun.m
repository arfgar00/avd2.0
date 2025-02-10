function = Ribfun(eqt,)

        I = obj.chord * (eqt)^3 / 12 + obj.chord * (eqt) * (obj.wingbox_depth/2)^2;
        
        function C = computeCrushLoad(obj)
            I = obj.secondMomentOfArea();
            C = (obj.bending_moment^2 * obj.rib_spacing * obj.wingbox_depth * ...
                (eqt) * obj.chord / 2) / ...
                (obj.young_mod * 10^9 * I^2);
        end
        
        function CS = computeCrushStress(obj, design_rib_thickness)
            crush_load = obj.computeCrushLoad();
            CS = (crush_load / (design_rib_thickness / 1000)) / obj.chord / 10^6; % MPa
        end
        
        function thickness_buckling = computeBucklingThickness(obj)
            crush_load = obj.computeCrushLoad();
            thickness_buckling = ((crush_load / obj.chord) / (3.62 * (obj.young_mod * 10^9)) * obj.wingbox_depth^2)^(1/3) * 1000;
        end
        
        function BS = computeBucklingStress(obj)
            thickness_buckling = obj.computeBucklingThickness();
            crush_load = obj.computeCrushLoad();
            BS = (crush_load / (thickness_buckling / 1000)) / obj.chord / 10^6; % MPa
        end
        
        function thickness_yield = computeYieldThickness(obj)
            crush_load = obj.computeCrushLoad();
            thickness_yield = (crush_load / (obj.yield_stress * 10^6)) / obj.chord * 1000;
        end
end
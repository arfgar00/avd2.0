classdef WingboxAnalysis
    properties
        bending_moment
        chord
        equivalent_panel_thickness
        rib_spacing
        wingbox_depth
        young_mod
        yield_stress = 350; % MPa
    end
    
    methods
        function obj = WingboxAnalysis(bm, c, ept, rs, wbd, ym)
            obj.bending_moment = bm;
            obj.chord = c;
            obj.equivalent_panel_thickness = ept;
            obj.rib_spacing = rs;
            obj.wingbox_depth = wbd;
            obj.young_mod = ym;
        end
        
        function I = secondMomentOfArea(obj)
            I = obj.chord * (obj.equivalent_panel_thickness/1000)^3 / 12 + ...
                obj.chord * (obj.equivalent_panel_thickness/1000) * (obj.wingbox_depth/2)^2;
        end
        
        function C = computeCrushLoad(obj)
            I = obj.secondMomentOfArea();
            C = (obj.bending_moment^2 * obj.rib_spacing * obj.wingbox_depth * ...
                (obj.equivalent_panel_thickness/1000) * obj.chord / 2) / ...
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
        
        function displayResults(obj)
            fprintf('Second Moment of Area: %.6e m^4\n', obj.secondMomentOfArea());
            fprintf('Crush Load: %.6f N\n', obj.computeCrushLoad());
            fprintf('Crush Stress (1mm rib thickness): %.6f MPa\n', obj.computeCrushStress(1));
            fprintf('Design Rib Thickness (Buckling): %.6f mm\n', obj.computeBucklingThickness());
            fprintf('Buckling Stress: %.6f MPa\n', obj.computeBucklingStress());
            fprintf('Design Rib Thickness (Yield): %.6f mm\n', obj.computeYieldThickness());
        end
    end
end

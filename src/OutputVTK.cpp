#include "OutputVTK.h"
#include "Forest.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"
#include "vtkHexahedron.h"

namespace trinurbs
{
    void OutputVTK::outputGeometry(const Forest& f) const
    {
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1; // number of cells in each parametric direction
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        
        uint sample_offset = 0;
        
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto e = f.bezierElement(i);
            
            uint count = 0;
            for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                const ParamPt samplept = isamplept.getCurrentPt();
                const Point3D phys_coord = e->eval(samplept.s, samplept.t);
                points->InsertPoint(sample_offset + count, phys_coord.data());
                ++count;
            }
            for( uint t = 0; t < ncell; ++t ) {
                for( uint s = 0; s < ncell; ++s ) {
                    vtkSmartPointer< vtkCell > cell = vtkQuad::New();
                    cell->GetPointIds()->SetId(0, sample_offset + t * nsample + s );
                    cell->GetPointIds()->SetId(1, sample_offset + t * nsample + s + 1 );
                    cell->GetPointIds()->SetId(2, sample_offset + ( t + 1 ) * nsample + s + 1 );
                    cell->GetPointIds()->SetId(3, sample_offset + ( t + 1 ) * nsample + s );
                    grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds() );
                }
            }
            sample_offset += nsample * nsample;
        }
        grid->SetPoints(points);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_geometry.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        if(!writer->Write())
            error( "Cannot write vtk file" );
    }
    
    void OutputVTK::outputNodalField(const Forest& f,
                                     const std::string& fieldname,
                                     const std::vector<double>& soln) const
    {
        
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1; // number of cells in each parametric direction
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        vtkSmartPointer<vtkDoubleArray> solndata = vtkDoubleArray::New();
        solndata->SetNumberOfComponents(3);
        solndata->SetName(fieldname.c_str());
        solndata->SetComponentName(0, "real");
        solndata->SetComponentName(1, "imag");
        solndata->SetComponentName(2, "abs");
        
        uint sample_offset = 0;
        
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto e = f.bezierElement(i);
            
            const auto gbasisivec = e->globalBasisFuncI();
            uint count = 0;
            for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept) {
                const ParamPt samplept = isamplept.getCurrentPt();
                const Point3D phys_coord = e->eval(samplept.s, samplept.t);
                points->InsertPoint(sample_offset + count, phys_coord.data());
                
                // temp code
                //                const Point3D p(0.0, 1.0, 0.0);
                //                const Point3D k(5.0, 0.0, 0.0);
                //
                //                std::vector<std::complex<double>> result{p[0], p[1], p[2]};
                //                const auto wave = std::exp(-std::complex<double>(0.0, dot(k, phys_coord)));
                //                for(auto& r : result)
                //                    r *= wave;
                //
                //                const auto& exact_val = result;
                //
                //                const auto& t1 = e->tangent(samplept.s, samplept.t, nurbs::ParamDir::S);
                //                const auto& t2 = e->tangent(samplept.s, samplept.t, nurbs::ParamDir::T);
                //                const double jpiola = nurbs::cross(t1, t2).length();
                //
                //                DoubleVecVec j;
                //                j.push_back(t1.asVec());
                //                j.push_back(t2.asVec());
                //                j.push_back(
                //                            {
                //                                1.0/jpiola * (j[0][1] * j[1][2] - j[0][2] * j[1][1]),
                //                                1.0/jpiola * (j[0][2] * j[1][0] - j[0][0] * j[1][2]),
                //                                1.0/jpiola * (j[0][0] * j[1][1] - j[0][1] * j[1][0])
                //                            });
                //                auto jinv = inv3x3Mat(j);
                //                for(auto& row : jinv)
                //                    row.erase(row.begin() + 2);
                //
                //
                //                std::vector<std::complex<double>> fexact(2);
                //                std::vector<std::complex<double>> f_h(2);
                //                for(size_t i = 0; i < 2; ++i)
                //                    for(size_t j = 0; j < 3; ++j)
                //                        fexact[i] += jpiola * jinv[j][i] * exact_val[j];
                //
                //                const auto basisvec = e->localBasis(samplept.s, samplept.t);
                //
                //                std::complex<double> val;
                //                for(uint ibasis = 0; ibasis < basisvec.size(); ++ibasis)
                //                    val += soln[gbasisivec[ibasis]] * basisvec[ibasis];
                
                // end temp code
                
                // now interpolate solution
                const auto basisvec = e->basis(samplept.s, samplept.t);
                
                std::complex<double> val;
                for(uint ibasis = 0; ibasis < basisvec.size(); ++ibasis)
                    val += soln[gbasisivec[ibasis]] * basisvec[ibasis];
                
                
                solndata->InsertComponent(sample_offset + count, 0, val.real());
                solndata->InsertComponent(sample_offset + count, 1, val.imag());
                solndata->InsertComponent(sample_offset + count, 2, std::abs(val));
                ++count;
            }
            for( uint t = 0; t < ncell; ++t ) {
                for( uint s = 0; s < ncell; ++s ) {
                    vtkSmartPointer< vtkCell > cell = vtkQuad::New();
                    cell->GetPointIds()->SetId(0, sample_offset + t * nsample + s );
                    cell->GetPointIds()->SetId(1, sample_offset + t * nsample + s + 1 );
                    cell->GetPointIds()->SetId(2, sample_offset + ( t + 1 ) * nsample + s + 1 );
                    cell->GetPointIds()->SetId(3, sample_offset + ( t + 1 ) * nsample + s );
                    grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds() );
                }
            }
            sample_offset += nsample * nsample;
        }
        grid->SetPoints(points);
        grid->GetPointData()->AddArray(solndata);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_complex_soln.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        if(!writer->Write())
            error( "Cannot write vtk file" );
    }
    
}
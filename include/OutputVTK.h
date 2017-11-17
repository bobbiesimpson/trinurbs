#ifndef TRINURBS_OUTPUTVTK_H
#define TRINURBS_OUTPUTVTK_H

#include <string>
#include <vector>
#include <mutex>

#include "base.h"
#include "IParentSample.h"
#include "Point3D.h"
#include "Forest.h"
#include "MultiscaleForest.h"

#include "OutputVTK.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"
#include "vtkHexahedron.h"
#include "vtkCellArray.h"
#include "vtkDataSetMapper.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

namespace trinurbs
{
    /// Forward declarations
    class Forest;
    class MultiscaleForest;
    
    /// Class responsible for output to VTK format
    class OutputVTK {
        
    public:
        
        /// Default constructor
        OutputVTK() : OutputVTK("unnamed_output") {}
        
        /// Construct with filename
        OutputVTK(const std::string& f,
                  const uint nsample = DEFAULT_NGRID_PTS)
        :
            mFilename(f),
            mSamplePtN(nsample) {}
        
        /// output the geometry of a forest to VTK
        void outputForestGeometry(const Forest& f) const;
        
        /// output geometry of multiscale forest
        void outputMultiscaleForestGeometry(const MultiscaleForest& f) const;
        
        /// Write a complex nodal field to a vtu file
        /// if nlocaldof > 1 then we assume the soln is ordered
        /// as e.g. [u^1_x u^1_y u^1_z u^2_x u^2_y ...]
        void outputNodalField(const Forest& f,
                              const std::string& fieldname,
                              const std::vector<double>& soln,
                              const uint nlocaldof = 1) const;
        
        /// Output nodal field of a multiscale forest to vtu file
        template<typename S>
        void outputNodalField(const MultiscaleForest& f,
                              const std::string& fieldname,
                              const S& soln,
                              const uint nlocaldof = 1) const
        {
            const uint nsample = samplePtN();   // number of sample points in each parametric direction
            const uint ncell = nsample - 1;     // number of cells in each parametric direction
            
            // create vtk data structures
            vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
            vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
            
            vtkSmartPointer<vtkDoubleArray> solndata = vtkDoubleArray::New();
            
            solndata->SetNumberOfComponents(nlocaldof);
            solndata->SetName(fieldname.c_str());
            
            // set comonent names
            for(uint idof = 0; idof < nlocaldof; ++idof)
            {
                const std::string s = "component " + std::to_string(idof);
                solndata->SetComponentName(idof, s.c_str());
            }
            
            uint sample_offset = 0;
            
            // loop over element
            for(uint i = 0; i < f.elemN(); ++i)
            {
                const auto e = f.multiscaleBezierElement(i);
                const auto gbasisivec = e->globalBasisIVec();
                
                uint count = 0;
                
                // loop over sample points
                for(IParentSample isamplept(nsample); !isamplept.isDone(); ++isamplept)
                {
                    auto samplept = isamplept.getCurrentPt();
                    
                    // get physical coordinate of sample point
                    const Point3D phys_coord = e->eval(samplept.u, samplept.v, samplept.w);
                    points->InsertNextPoint(phys_coord[0], phys_coord[1], phys_coord[2]);
                    //points->InsertPoint(sample_offset + count, phys_coord.data());
                    
                    // now interpolate soluiton
                    const auto basisvec = e->basis(samplept.u,
                                                   samplept.v,
                                                   samplept.w);
                    
                    std::vector<double> val(nlocaldof, 0.0);
                    
                    for(uint ibasis = 0; ibasis < basisvec.size(); ++ibasis)
                        for(uint idof = 0; idof < nlocaldof; ++idof)
                            val[idof] += soln[gbasisivec[ibasis] * nlocaldof + idof] * basisvec[ibasis];
                    
                    for(uint idof = 0; idof < nlocaldof; ++idof)
                        solndata->InsertComponent(sample_offset + count, idof, val[idof]);
                    
                    ++count;
                }
                
                // now loop over cells and construct connectivity
                for(uint iw = 0; iw < ncell; ++iw)
                {
                    for(uint iv = 0; iv < ncell; ++iv)
                    {
                        for(uint iu = 0; iu < ncell; ++iu)
                        {
                            // create a hex
                            vtkSmartPointer<vtkHexahedron> cell = vtkSmartPointer<vtkHexahedron>::New();
                            
                            const uint index0 = iw * nsample * nsample + iv * nsample + iu;
                            const uint index2 = iw * nsample * nsample + (iv + 1) * nsample + iu + 1;
                            
                            const uint index4 = (iw + 1) * nsample * nsample + iv * nsample + iu;
                            const uint index6 = (iw + 1) * nsample * nsample + (iv + 1) * nsample + iu + 1;
                            
                            cell->GetPointIds()->SetId(0, sample_offset + index0);
                            cell->GetPointIds()->SetId(1, sample_offset + index0 + 1);
                            cell->GetPointIds()->SetId(2, sample_offset + index2);
                            cell->GetPointIds()->SetId(3, sample_offset + index2 - 1);
                            cell->GetPointIds()->SetId(4, sample_offset + index4);
                            cell->GetPointIds()->SetId(5, sample_offset + index4 + 1 );
                            cell->GetPointIds()->SetId(6, sample_offset + index6);
                            cell->GetPointIds()->SetId(7, sample_offset + index6 - 1);
                            
                            cellarray->InsertNextCell(cell);
                        }
                    }
                }
                sample_offset += nsample * nsample * nsample;
            }
            grid->SetPoints(points);
            grid->SetCells(VTK_HEXAHEDRON, cellarray);
            grid->GetPointData()->AddArray(solndata);
            
            vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
            
            const std::string fname = filename() + "_nodal_soln.vtu";
            
            writer->SetFileName(fname.c_str());
            writer->SetInputData(grid);
            
            if(!writer->Write())
                error( "Cannot write vtk file" );
        }
        
        /// Sample point number setter
        void setSamplePtN(const uint n)
        {
            if(n < 1)
                error("Must specifiy non-zero sample points for output");
            mSamplePtN = n;
        }
        
        /// Filename getter
        const std::string& filename() const
        {
            return mFilename;
        }
        
        /// Filename setter
        void setFilename(const std::string& newfile)
        {
            mFilename = newfile;
        }
        
        /// Sample point number getter
        uint samplePtN() const { return mSamplePtN; }
        
    private:
        
        /// Filename we write to
        std::string mFilename;
        
        /// Number of sample points
        uint mSamplePtN;
        
        /// Get the mutex
        std::mutex& mutex() const
        {
            return mMutex;
        }
        
        /// Mutex for this class
        mutable std::mutex mMutex;
    };
    
}

#endif
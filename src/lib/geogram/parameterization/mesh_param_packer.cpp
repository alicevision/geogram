/*
 *  Copyright (c) 2012-2016, Inria, the ALICE project.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram/parameterization/mesh_param_packer.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/logger.h>
#include <geogram/numerics/matrix_util.h>
#include <algorithm>
#include <math.h>

namespace {
    using namespace GEO;

    class AverageDirection2d {
    public:
	
        AverageDirection2d()  {
	    result_is_valid_ = false;
	}
	
        void begin() {
	    for(int i=0; i<3; i++) {
		M_[i] = 0.0 ;
	    }
	    result_is_valid_ = false;	    
	}
	
        void add_vector(const vec2& e) {
	    M_[0] += e.x * e.x ;
	    M_[1] += e.x * e.y ;
	    M_[2] += e.y * e.y ;
	}
	
        void end() {
	    double eigen_vectors[4] ;
	    double eigen_values[2] ;
	    MatrixUtil::semi_definite_symmetric_eigen(
		M_, 2, eigen_vectors, eigen_values
	    );
	    int k = 0 ;
	    result_ = vec2(
		eigen_vectors[2*k],
		eigen_vectors[2*k+1]
	    );
	    result_is_valid_ = true ;
	}
	
        const vec2& average_direction() const { 
            geo_assert(result_is_valid_) ;
            return result_ ; 
        }
	
    private:
        double M_[3] ;
        vec2 result_ ;
        bool result_is_valid_ ;
    } ;
    
    
}


namespace GEO {

    /**
     * \brief Stores a reference to a chart, its bounding box
     *  in parameter space, and its upper and lower horizons.
     */
    class ChartBBox {
    public:

	/**
	 * \brief ChartBBox constructor.
	 * \param[in] chart a pointer to the chart.
	 * \param[in] min , max lower-left and upper-right
	 *  corners of the parameter-space bounding box of
	 *  the chart.
	 */
        ChartBBox(
            Chart* chart, const vec2& min, const vec2& max
        ) : 
            min_(min), max_(max), chart_(chart),
	    min_func_(nullptr), max_func_(nullptr),
            nb_steps_(0) {
        }

	/**
	 * \brief ChartBBox destructor.
	 */
        ~ChartBBox(){
            free();
            chart_ = nullptr;
        }

	/**
	 * \brief ChartBBox copy-constructor.
	 * \param[in] rhs a const reference to the ChartBBox to be copied.
	 */
        ChartBBox(const ChartBBox& rhs) {
            copy(rhs);
        }

	/**
	 * \brief ChartBBox affectation.
	 * \param[in] rhs a const reference to the ChartBBox to be copied.
	 * \return A reference to this ChartBBox after affectation.
	 */
        ChartBBox& operator=(const ChartBBox& rhs) {
            if(&rhs != this) {
                free();
                copy(rhs);
            }
            return *this;    
        }

	/**
	 * \brief Initializes the upper and lower horizons.
	 * \param[in] step pixel-size in parameter space
	 * \param[in] margin margin size in parameter space
	 * \param[in] margin_width_in_pixels margin size in pixels
	 * \param[in] tex_coord a 2d vector attribute attached to 
	 *  the facet corners of the mesh with the texture coordinates.
	 * \param[in] chart_attr an attribute attached to the facets of
	 *  the mesh with the id of the chart the facet belongs to.
	 */
        void init_max_and_min_func(
            double step, double margin, index_t margin_width_in_pixels,
	    Attribute<double>& tex_coord, Attribute<index_t>& chart_attr
        ) {

            if(min_func_ != nullptr) {
                delete[] min_func_;
            }
            if(max_func_ != nullptr) {
                delete[] max_func_;
            }

            nb_steps_ = int(1.0 + (max_.x - min_.x) / step );
            min_func_ = new double[nb_steps_];
            max_func_ = new double[nb_steps_];
            for (int i=0;i<nb_steps_;i++){
                min_func(i) = max_.y - min_.y;
                max_func(i) = 0;
            }

	    for(index_t ff=0; ff<chart_->facets.size(); ++ff) {
		index_t f = chart_->facets[ff];
		for(
		    index_t c1 = chart_->mesh.facets.corners_begin(f);
		    c1 < chart_->mesh.facets.corners_end(f); ++c1
		) {
		    index_t neighf = chart_->mesh.facet_corners.adjacent_facet(c1);
		    
		    if(neighf != NO_FACET && chart_attr[neighf] == chart_->id) {
			continue;
		    }
		    
		    index_t c2 = chart_->mesh.facets.next_corner_around_facet(f,c1);
		    vec2 p1(tex_coord[2*c1], tex_coord[2*c1+1]);
		    vec2 p2(tex_coord[2*c2], tex_coord[2*c2+1]);

		    if(p2.x == p1.x) {
			continue;
		    }
		    
		    if(p2.x < p1.x) {
			std::swap(p1,p2);
		    }

		    p1.x -= margin_width_in_pixels * step;
		    p2.x += margin_width_in_pixels * step;
		    
		    double a = (p2.y - p1.y) / (p2.x - p1.x);
		    for(double x = p1.x; x<p2.x; x+=step) {
			double y = p1.y + a * (x - p1.x);
			int cX = int((x - min_.x) / step);
                        if (cX>=0 && cX<nb_steps_) {
                            min_func(cX) = std::min(min_func(cX), y - margin);
                            max_func(cX) = std::max(max_func(cX), y + margin);
			}
		    }
		}
	    }
        }

        double& max_func(int i) {
            geo_assert(i>=0 && i<nb_steps_);
            return max_func_[i];
        }

        double& min_func(int i) {
            geo_assert(i>=0 && i<nb_steps_);
            return min_func_[i];
        }

        vec2 size() const {
            return max_-min_;
        }

	/**
	 * \brief Gets the area of the bounding box.
	 * \return the area of the bounding box in parameter space.
	 */
        double area() const {
            vec2 s = size();
            return s.x*s.y;
        }

	/**
	 * \brief Applies a translation vector to the texture coordinates.
	 * \param[in] v the translation vector.
	 * \param[in,out] tex_coord a dimension 2 vector attribute attached to
	 *  the facet corners of the mesh with the texture coordinates.
	 */
        void translate(const vec2& v, Attribute<double>& tex_coord){
            min_=min_+v;
            max_=max_+v; 
	    for(index_t ff=0; ff<chart_->facets.size(); ++ff) {
		index_t f = chart_->facets[ff];
		for(
		    index_t c = chart_->mesh.facets.corners_begin(f);
		    c < chart_->mesh.facets.corners_end(f); ++c
		) {
		    tex_coord[2*c] += v.x;
		    tex_coord[2*c+1] += v.y;
		}
	    }
        }

	/**
	 * \brief Gets the lower-left corner of the bounding box.
	 * \return A const reference to the parametric-space 
	 *  coordinates of the lower-left corner of the bounding box.
	 */
        const vec2& min() const {
	    return min_;
	}

	/**
	 * \brief Gets the upper-right corner of the bounding box.
	 * \return A const reference to the parametric-space 
	 *  coordinates of the upper-right corner of the bounding box.
	 */
	const vec2& max() const {
	    return max_;
	}

	/**
	 * \brief Gets the lower-left corner of the bounding box.
	 * \return A modifiable reference to the parametric-space 
	 *  coordinates of the lower-left corner of the bounding box.
	 */
        vec2& min() {
	    return min_;
	}

	/**
	 * \brief Gets the upper-right corner of the bounding box.
	 * \return A modifiable reference to the parametric-space 
	 *  coordinates of the upper-right corner of the bounding box.
	 */
	vec2& max() {
	    return max_;
	}

	/**
	 * \brief Gets the chart.
	 * \return a const pointer to the chart.
	 */
        const Chart* chart() const {
	    return chart_;
	}


    protected:

	/**
	 * \brief Copies a ChartBBox.
	 * \param[in] rhs a const reference to the ChartBBox 
	 *  to be copied.
	 */
        void copy(const ChartBBox& rhs) {
            chart_ = rhs.chart_;
            nb_steps_ = rhs.nb_steps_;

            if(rhs.min_func_ != nullptr) {
                min_func_ = new double[nb_steps_];
                max_func_ = new double[nb_steps_];
                for(int i=0; i<nb_steps_; i++) {
                    min_func_[i] = rhs.min_func_[i];
                    max_func_[i] = rhs.max_func_[i];
                }
            } else {
                geo_assert(rhs.max_func_ == nullptr);
                min_func_ = nullptr;
                max_func_ = nullptr;
            }
            min_ = rhs.min_;
            max_ = rhs.max_;
        }

	/**
	 * \brief Releases the memory allocated by
	 *  this ChartBBox.
	 */
        void free() {
            delete[] min_func_;
            delete[] max_func_;
            min_func_ = nullptr;
            max_func_ = nullptr;
            nb_steps_ = 0;
        }

    private: 
        vec2 min_, max_;
        Chart* chart_;
        double* min_func_;
        double* max_func_;
        int nb_steps_;
    };

    /***********************************************************/

    /**
     * \brief Internal implementation of the packing algorithm.
     * \details By Nicolas Ray. The algorithm is inspired by 
     *  the Tetris game, and iteratively places the charts from 
     *  bottom to top.
     */
    class TetrisPacker {
    public :

	/**
	 * \brief TetrisPacker constructor.
	 * \param[in] mesh a reference the surface mesh to be packed.
	 *  It needs to have a vector attribute of dimension 2 attached
	 *  to the facet corners and named "tex_coord", as well as a 
	 *  facet attribute named "chart".
	 */
        TetrisPacker(Mesh& mesh) {
            nb_xpos_ = 1024;
            height_ = new double[nb_xpos_];
            image_size_in_pixels_  = 1024;
            margin_width_in_pixels_ = 4;
	    tex_coord_.bind_if_is_defined(
		mesh.facet_corners.attributes(), "tex_coord"
	    );
	    geo_assert(tex_coord_.is_bound() && tex_coord_.dimension() == 2);
	    chart_attr_.bind_if_is_defined(
		mesh.facets.attributes(), "chart"
	    );
	    geo_assert(chart_attr_.is_bound());
        }

	/**
	 * \brief TetrisPacker destructor.
	 */
        ~TetrisPacker() {
            delete[] height_;
            height_ = nullptr;
        }

        void set_image_size_in_pixels(index_t size) {
            image_size_in_pixels_ = size;
        }

        index_t margin_width_in_pixels() const {
	    return margin_width_in_pixels_;
	}

        void set_margin_width_in_pixels(index_t width) {
            margin_width_in_pixels_ = width;
        } 
        
        void add(const ChartBBox& r){
            data_.push_back(r);
        }

        /**
	 * \brief Compares two ChartBBox objects.
	 * \details Used by the packing algorithm to sort
	 *  the boxes.
	 * \param[in] b0 , b1 const references to the
	 *  two ChartBBox objects to be compared.
	 * \retval true if \p b0 is taller than \p b1.
	 * \retval false otherwise.
	 */
        static bool compare(
            const ChartBBox& b0, const ChartBBox& b1
        ) {
            return (b0.size().y > b1.size().y);
        }
        
        void add_margin(double margin_size) {
            for (unsigned int i=0;i<data_.size();i++){
                data_[i].translate(
                    vec2(
                        -data_[i].min().x + margin_size,
                        -data_[i].min().y + margin_size
		    ),
		    tex_coord_
                );
                data_[i].max() = data_[i].max() + 
                    vec2( 2.0 * margin_size, 2.0 * margin_size);
                geo_assert(data_[i].max().x>0);
                geo_assert(data_[i].max().y>0);
                data_[i].min() = vec2(0,0);
            }
        }
        
        double max_height() {
            double result = 0;
            for (unsigned int i=0;i<data_.size();i++) {
                result = std::max(result, data_[i].max().y);
            }
            return result;
        }

        void recursive_apply() {
            // compute margin
            double area = 0;      
            for (unsigned int numrect = 0; numrect <data_.size();numrect++ ) {
                area += data_[numrect].area();
            }
            double margin = 
                (::sqrt(area) / image_size_in_pixels_) * 
                margin_width_in_pixels_;
            add_margin(margin); 

            // find a first solution
            apply(margin);
            double scoreSup = max_height();
            double borneSup = width_;
            double decalborne = 0.5 * ::sqrt(scoreSup * width_);

            // dichotomy
            for  (int i=0;i<10;i++ ){
                double new_borne = borneSup - decalborne;
                apply(margin,new_borne);
                double max = max_height();
                if (max < new_borne){
                    borneSup = borneSup - decalborne;
                }
                decalborne/=2;
            }
            apply(margin, borneSup);
        }


        void apply(double margin, double width =-1) {
            width_ = width;
            for (unsigned int i=0; i<data_.size(); i++) {
                data_[i].translate(
                    vec2(
                        -data_[i].min().x,
                        -data_[i].min().y
		    ),
		    tex_coord_
                );
            }

            // sort ChartBBoxes by their heights ...
            std::sort(data_.begin(),data_.end(),compare);            
            
            { //find the best width
                double max_bbox_width = 0;
                for (unsigned int i=0; i<data_.size(); i++) {
                    max_bbox_width= std::max(
                        max_bbox_width, data_[i].size().x
                    );
                }
                
                if ( width == -1){
                    // try to have a square : width_ = sqrt(area);
                    double area =0;
                
                    for (unsigned int i=0; i<data_.size(); i++){
                        area += data_[i].area();
                    }

                    width_ = ::sqrt(area) *1.1;
                

                    // resize if a square seems to be too bad 
                    // (the first piece is too high)

                    if (data_[0].size().y > width_) {
                        width_ = area / data_[0].size().y;
                    }
                }

                // be sure all surface can fit in the width
                width_ = std::max(
                    max_bbox_width * (double(nb_xpos_ + 2) / double(nb_xpos_)),
                    width_
                );	
            }
            // set the step depending on the width and the discretisation
            step_ = width_ / nb_xpos_;


            // init local min and max height functions
	    for (unsigned int numrect = 0; numrect <data_.size(); numrect++) {
                data_[numrect].init_max_and_min_func(
                    step_, margin, margin_width_in_pixels_, tex_coord_, chart_attr_
                );
            }

            // init global height function
            for (int x=0; x<nb_xpos_; x++) {
                height(x)=0;
            }
            
            for (unsigned int i=0;i<data_.size();i++) {
                place(data_[i]);
            }
        }
        
        
    private:

        double& height(int x) {
            geo_assert(x >= 0 && x < nb_xpos_);
            return height_[x];
        }

        const double& height(int x) const {
            geo_assert(x >= 0 && x < nb_xpos_);
            return height_[x];
        }
        

	/**
	 * \brief Inserts a new chart.
	 * \details Follows the "tetris" strategy.
	 * \param[in] rect the chart to be inserted.
	 */
        void place(ChartBBox& rect) {

            const int width_in_pas = int ( (rect.size().x / step_) + 1 );
            
            // find the best position
            int bestXPos=0;
            double bestHeight = Numeric::max_float64();
            double bestFreeArea = Numeric::max_float64();
            
            for (int x=0;x<nb_xpos_ - width_in_pas;x++) {
                
                double localHeight = max_height(x,width_in_pas,rect);
                double localFreeArea = 
                    free_area(x,width_in_pas,localHeight)
                    + localHeight * rect.size().x;
                
                // the best position criterion is the position that 
                // minimize area lost
                // area could be lost under and upper the bbox
                if (localFreeArea < bestFreeArea){
                    bestXPos = x;
                    bestHeight = localHeight;
                    bestFreeArea = localFreeArea;
                }
            }

            // be sure we have a solution
            geo_assert(bestHeight != Numeric::max_float64());
            
            // place the rectangle
            rect.translate(vec2(bestXPos*step_,bestHeight),tex_coord_);
            for (int i=bestXPos;i<bestXPos + width_in_pas;i++){
                height(i) =  rect.min().y + rect.max_func(i-bestXPos);
            }
        }
        
        // find the max height in the range  [xpos,  xpos + width]
        double max_height(int xpos, int width, ChartBBox& rect){
            double result = 0;
            for (int i=xpos; i<xpos+width; i++) {
                result = std::max( result, height(i)-rect.min_func(i-xpos) );
            }
            return result;
        }
        
        // find the free area in the range under height_max in range
        // [xpos,  xpos + width]

        double free_area(int xpos, int width, double height_max) {
            double result =0;
            for (int i=xpos; i<xpos+width; i++) {
                result += step_ * (height_max - height(i));
            }
            return result;
        }

    private:
        index_t image_size_in_pixels_;
        index_t margin_width_in_pixels_;

        int nb_xpos_;
        double width_;
        double step_;
        double *height_;

        std::vector<ChartBBox> data_;

	Attribute<double> tex_coord_;
	Attribute<index_t> chart_attr_;
    };


/***********************************************************************/


    Packer::Packer() {
        image_size_in_pixels_ = 1024;
        margin_width_in_pixels_ = 4;
    }
  

    // There where some problems with nan (Not a number),
    // this code is used to track such problems (seems to
    // be ok now)
    static bool chart_is_ok(Chart& chart, Attribute<double>& tex_coord) {
	for(index_t ff=0; ff<chart.facets.size(); ++ff) {
	    index_t f = chart.facets[ff];
	    for(
		index_t c = chart.mesh.facets.corners_begin(f);
		c < chart.mesh.facets.corners_end(f); ++c
	    ) {
		if(Numeric::is_nan(tex_coord[2*c])) {
		    return false;
		} 
		if(Numeric::is_nan(tex_coord[2*c+1])) {
		    return false;
		} 
	    }
	}
	return true;
    }


    void Packer::pack_surface(Mesh& mesh) {
	tex_coord_.bind_if_is_defined(
	    mesh.facet_corners.attributes(), "tex_coord"
	);
	geo_assert(tex_coord_.is_bound() && tex_coord_.dimension() == 2);
	chart_attr_.bind_if_is_defined(
	    mesh.facets.attributes(), "chart"
	);
	geo_assert(chart_attr_.is_bound());

	// Get the charts
	vector<Chart> charts;
	index_t nb_charts=0;
	for(index_t f=0; f<mesh.facets.nb(); ++f) {
	    nb_charts = std::max(nb_charts, chart_attr_[f]);
	}
	++nb_charts;

	for(index_t i=0; i<nb_charts; ++i) {
	    charts.push_back(Chart(mesh, i));
	}
	
	for(index_t f=0; f<mesh.facets.nb(); ++f) {
	    charts[chart_attr_[f]].facets.push_back(f);
	}

	Logger::out("Packer") << "Packing " << charts.size() << " charts" << std::endl;
	
        // Sanity check
	for(index_t i=0; i<nb_charts; ++i) {
            if(!chart_is_ok(charts[i], tex_coord_)) {
		for(index_t ff=0; ff<charts[i].facets.size(); ++ff) {
		    index_t f = charts[i].facets[ff];
		    for(
			index_t c = mesh.facets.corners_begin(f);
			c<mesh.facets.corners_end(f); ++c
		    ) {
			tex_coord_[2*c] = 0.0;
			tex_coord_[2*c+1] = 0.0;
		    }
		}
            }
        }

	pack_charts(charts);

	// Normalize tex coords in [0,1] x [0,1]
	{
	    double u_min, v_min, u_max, v_max;
	    Geom::get_mesh_bbox_2d(
		mesh, tex_coord_, u_min, v_min, u_max, v_max
	    );
	    double l = std::max(u_max - u_min, v_max - v_min);
	    if(l > 1e-6) {
		for(index_t c=0; c<mesh.facet_corners.nb(); ++c) {
		    tex_coord_[2*c]   = (tex_coord_[2*c] - u_min) / l;
		    tex_coord_[2*c+1] = (tex_coord_[2*c+1] - v_min) / l;		    
		}
	    }
	}

	tex_coord_.unbind();
	chart_attr_.unbind();
    }


    void Packer::pack_charts(vector<Chart>& charts) {
        
        Logger::out("Packer") 
            << "nb components:" << charts.size() << std::endl;  

	if(charts.size() == 0) {
	    return;
	}
	
        Mesh& mesh = charts[0].mesh;
        total_area_3d_ = Geom::mesh_area(mesh);

        for(index_t i=0; i<charts.size(); ++i) {
            geo_assert(chart_is_ok(charts[i], tex_coord_));
            normalize_chart(charts[i]);
            geo_assert(chart_is_ok(charts[i], tex_coord_));
        }

        // use the tetris packer (more efficient for large dataset)
        // set some application dependent const
        TetrisPacker pack(mesh);  
        pack.set_image_size_in_pixels(image_size_in_pixels());
        pack.set_margin_width_in_pixels(margin_width_in_pixels());
        double area = 0;
        for(index_t i=0; i<charts.size(); ++i) {
            geo_assert(chart_is_ok(charts[i], tex_coord_));
            double u_min, v_min, u_max, v_max;
	    Geom::get_chart_bbox_2d(charts[i], tex_coord_, u_min, v_min, u_max, v_max);
            
            geo_assert(!Numeric::is_nan(u_min));
            geo_assert(!Numeric::is_nan(v_min));
            geo_assert(!Numeric::is_nan(u_max));
            geo_assert(!Numeric::is_nan(v_max));
            geo_assert(u_max >= u_min);
            geo_assert(v_max >= v_min);
            
            area += (v_max - v_min) * (u_max - u_min);
            
            ChartBBox r(
                &charts[i],
                vec2(u_min,v_min) ,
                vec2(u_max,v_max)
            );
            pack.add(r);
        }
        
        
        pack.recursive_apply();

	{
	    double total_area = Geom::mesh_area_2d(mesh, tex_coord_);
	    double u_min, v_min, u_max, v_max;
	    Geom::get_mesh_bbox_2d(mesh, tex_coord_, u_min, v_min, u_max, v_max);
            double bbox_area = (u_max - u_min)*(v_max-v_min);
            double filling_ratio = total_area / bbox_area;
            Logger::out("Packer") << "BBox area:"  << bbox_area << std::endl;
            Logger::out("Packer") << "Filling ratio:" 
                                  << filling_ratio << std::endl;
        }
    }

    void Packer::normalize_chart(Chart& chart) {

	// TODO: choose best axes, rotate chart.
	vec2 U(1.0, 0.0);
	vec2 V(0.0, 1.0);
	AverageDirection2d dir;
	dir.begin();
	for(index_t ff=0; ff<chart.facets.size(); ++ff) {
	    index_t f = chart.facets[ff];
	    for(index_t c1=chart.mesh.facets.corners_begin(f);
		c1 < chart.mesh.facets.corners_end(f); ++c1) {
		index_t adj_f = chart.mesh.facet_corners.adjacent_facet(c1);
		if(adj_f == NO_FACET || chart_attr_[adj_f] != chart.id) {
		    index_t c2 = chart.mesh.facets.next_corner_around_facet(f,c1);
		    vec2 uv1(tex_coord_[2*c1], tex_coord_[2*c1+1]);
		    vec2 uv2(tex_coord_[2*c2], tex_coord_[2*c2+1]);
		    dir.add_vector(uv2-uv1);
		}
	    }
	}
	dir.end();
	vec2 W=dir.average_direction();
	if(!Geom::has_nan(W)) {
	    W = normalize(W);
	    V = W;
	    U = vec2(-V.y, V.x);
	}
	
	for(index_t ff=0; ff<chart.facets.size(); ++ff) {
	    index_t f = chart.facets[ff];
	    for(index_t c=chart.mesh.facets.corners_begin(f);
		c < chart.mesh.facets.corners_end(f); ++c) {
		vec2 uv(tex_coord_[2*c], tex_coord_[2*c+1]);
		tex_coord_[2*c] = dot(uv,U);
		tex_coord_[2*c+1] = dot(uv,V);
	    }
	}
	
        double area3d = Geom::chart_area(chart);
        double area2d = Geom::chart_area_2d(chart, tex_coord_);
        double factor = 1.0;
        if(::fabs(area2d) > 1e-30) {
            factor = ::sqrt(area3d/area2d);
        } else {
            factor = 0.0;
        }

        geo_assert(!Numeric::is_nan(area2d));
        geo_assert(!Numeric::is_nan(area3d));
        geo_assert(!Numeric::is_nan(factor));

	for(index_t ff=0; ff<chart.facets.size(); ++ff) {
	    index_t f = chart.facets[ff];
	    for(index_t c=chart.mesh.facets.corners_begin(f);
		c < chart.mesh.facets.corners_end(f); ++c) {
		tex_coord_[2*c] *= factor;
		tex_coord_[2*c+1] *= factor;
	    }
	}
    }
}


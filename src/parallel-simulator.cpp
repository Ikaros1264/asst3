#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>
#include <omp.h>
// TASK 2

// NOTE: You may modify this class definition as you see fit, as long as the
// class name, and type of simulateStep and buildAccelerationStructure remain
// the same. You may modify any code outside this class unless otherwise
// specified.

const int QuadTreeLeafSize = 8;
class ParallelNBodySimulator : public INBodySimulator {
public:
  // TODO: implement a function that builds and returns a quadtree containing
  // particles. You do not have to preserve this function type.
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {
    std::unique_ptr<QuadTreeNode> This(new QuadTreeNode);
    This->isLeaf = 0;
    This->children[0] = nullptr; This->children[1] = nullptr; This->children[2] = nullptr; This->children[3] = nullptr; //set node to null

    if(particles.size() > QuadTreeLeafSize){
      Vec2 bmid((bmin.x + bmax.x) / 2, (bmin.y + bmax.y) / 2);
      std::vector<Particle> part[4];

      for(auto p : particles) {
        if(p.position.x>=bmin.x && p.position.x<=bmid.x && p.position.y>=bmin.y && p.position.y<=bmid.y){
          part[0].push_back(p);
        }
        else if(p.position.x>=bmid.x && p.position.x<=bmax.x && p.position.y>=bmin.y && p.position.y<=bmid.y){
          part[1].push_back(p);
        }
        else if(p.position.x>=bmin.x && p.position.x<=bmid.x && p.position.y>=bmid.y && p.position.y<=bmax.y){
          part[2].push_back(p);
        }
        else if(p.position.x>=bmid.x && p.position.x<=bmax.x && p.position.y>=bmid.y && p.position.y<=bmax.y){
          part[3].push_back(p);
        }
      }

      #pragma omp parallel sections
      {
        #pragma omp section
        This->children[0] = buildQuadTree(part[0],bmin,bmid);
        #pragma omp section
        This->children[1] = buildQuadTree(part[1],Vec2(bmid.x,bmin.y),Vec2(bmax.x,bmid.y));
        #pragma omp section
        This->children[2] = buildQuadTree(part[2],Vec2(bmin.x,bmid.y),Vec2(bmid.x,bmax.y));
        #pragma omp section
        This->children[3] = buildQuadTree(part[3],bmid,bmax);
      }

    }else{
      This->particles = particles;
      This->isLeaf = 1;
    }
    return This;
  }

  // Do not modify this function type.
  virtual std::unique_ptr<AccelerationStructure>
  buildAccelerationStructure(std::vector<Particle> &particles) {
    // build quad-tree
    auto quadTree = std::make_unique<QuadTree>();

    // find bounds
    Vec2 bmin(1e30f, 1e30f);
    Vec2 bmax(-1e30f, -1e30f);

    for (auto &p : particles) {
      bmin.x = fminf(bmin.x, p.position.x);
      bmin.y = fminf(bmin.y, p.position.y);
      bmax.x = fmaxf(bmax.x, p.position.x);
      bmax.y = fmaxf(bmax.y, p.position.y);
    }

    quadTree->bmin = bmin;
    quadTree->bmax = bmax;

    // build nodes
    quadTree->root = buildQuadTree(particles, bmin, bmax);
    if (!quadTree->checkTree()) {
      std::cout << "Your Tree has Error!" << std::endl;
    }

    return quadTree;
  }

  // Do not modify this function type.
  virtual void simulateStep(AccelerationStructure *accel,
                            std::vector<Particle> &particles,
                            std::vector<Particle> &newParticles,
                            StepParameters params) override {
    // TODO: implement parallel version of quad-tree accelerated n-body
    // simulation here, using quadTree as acceleration structure
    #pragma omp parallel for
    for(int i = 0; i<particles.size(); i++){
      Vec2 force = Vec2(0.0f, 0.0f);
      std::vector<Particle> nearby;
      accel->getParticles(nearby, particles[i].position, params.cullRadius);
      for(auto &n : nearby){
        force += computeForce(particles[i], n, params.cullRadius);
      }
      newParticles[i]=updateParticle(particles[i], force, params.deltaTime);
    }
  }
};

// Do not modify this function type.
std::unique_ptr<INBodySimulator> createParallelNBodySimulator() {
  return std::make_unique<ParallelNBodySimulator>();
}

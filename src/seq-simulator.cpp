#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>

// TASK 1

// NOTE: You may modify any of the contents of this file, but preserve all
// function types and names. You may add new functions if you believe they will
// be helpful.

const int QuadTreeLeafSize = 8;
class SequentialNBodySimulator : public INBodySimulator {
public:
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {
    // TODO: implement a function that builds and returns a quadtree containing
    // particles.
    std::unique_ptr<QuadTreeNode> This(new QuadTreeNode);
    This->isLeaf = 0;
    This->children[0] = nullptr; This->children[1] = nullptr; This->children[2] = nullptr; This->children[3] = nullptr; //set node to null
    for(int i = 0; i < particles.size(); i++) {
      auto p = particles[i];
      if(p.position.x>=bmin.x && p.position.x<=bmax.x && p.position.y>=bmin.y && p.position.y<=bmax.y){
        This->particles.push_back(p);
        std::swap(particles[i],particles.back());
        particles.pop_back();
        i--;
      }
    }
    if(This->particles.size() > QuadTreeLeafSize){
      Vec2 bmid((bmin.x + bmax.x) / 2, (bmin.y + bmax.y) / 2);
      //#pragma omp parallel sections
      {
        //#pragma omp section
        This->children[0] = buildQuadTree(This->particles,bmin,bmid);
        //#pragma omp section
        This->children[1] = buildQuadTree(This->particles,Vec2(bmid.x,bmin.y),Vec2(bmax.x,bmid.y));
        //#pragma omp section
        This->children[2] = buildQuadTree(This->particles,Vec2(bmin.x,bmid.y),Vec2(bmid.x,bmax.y));
        //#pragma omp section
        This->children[3] = buildQuadTree(This->particles,bmid,bmax);
      }
    }else{
      This->isLeaf = 1;
    }
    return This;
  }
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
    std::vector<Particle> cpy(particles);
    quadTree->root = buildQuadTree(particles, bmin, bmax);
    particles.swap(cpy);

    if (!quadTree->checkTree()) {
      std::cout << "Your Tree has Error!" << std::endl;
    }

    return quadTree;
  }
  virtual void simulateStep(AccelerationStructure *accel,
                            std::vector<Particle> &particles,
                            std::vector<Particle> &newParticles,
                            StepParameters params) override {
    // TODO: implement sequential version of quad-tree accelerated n-body
    // simulation here, using quadTree as acceleration structure
    newParticles.clear();
    for(auto &p : particles){
      Vec2 force = Vec2(0.0f, 0.0f);
      std::vector<Particle> nearby;
      accel->getParticles(nearby, p.position, params.cullRadius);
      for(auto &n : nearby){
        force += computeForce(p, n, params.cullRadius);
      }
      newParticles.push_back(updateParticle(p, force, params.deltaTime));
    }
    //newParticles.assign(particles.begin(),particles.end());
  }
};

std::unique_ptr<INBodySimulator> createSequentialNBodySimulator() {
  return std::make_unique<SequentialNBodySimulator>();
}

#ifndef EPIWORLDR_MODEL_BUILDER_H
#define EPIWORLDR_MODEL_BUILDER_H

#include "epiworld-common.h"

// One configurable model type for R recipes. Keeping the initializer factory
// on the native model preserves initial_states() after cloning and in threads.
class ModelBuilder : public epiworld::Model<> {
public:
  using Initializer = std::function<void(epiworld::Model<>*)>;
  using InitializerFactory = std::function<Initializer(std::vector<double>)>;
  InitializerFactory initializer_factory;

  epiworld::Model<> & initial_states(
      std::vector<double> proportions, std::vector<int> = {}
  ) override {
    if (initializer_factory)
      initial_states_fun = initializer_factory(std::move(proportions));
    return *this;
  }

  std::unique_ptr<epiworld::Model<>> clone_ptr() override {
    return std::make_unique<ModelBuilder>(*this);
  }
};

#endif

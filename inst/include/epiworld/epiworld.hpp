#include <vector>
#include <functional>
#include <memory>
#include <stdexcept>
#include <random>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <chrono>
#include <climits>
#include <cstdint>
#include <algorithm>
#include <regex>
#include <sstream>
#include <iomanip>
#include <set>
#include <type_traits>
#include <cassert>

#ifndef EPIWORLD_HPP
#define EPIWORLD_HPP

/* Versioning */
#define EPIWORLD_VERSION_MAJOR 0
#define EPIWORLD_VERSION_MINOR 8
#define EPIWORLD_VERSION_PATCH 2

static const int epiworld_version_major = EPIWORLD_VERSION_MAJOR;
static const int epiworld_version_minor = EPIWORLD_VERSION_MINOR;
static const int epiworld_version_patch = EPIWORLD_VERSION_PATCH;

namespace epiworld {

    #include "config.hpp"
    #include "epiworld-macros.hpp"

    #include "misc.hpp"
    #include "progress.hpp"

    #include "modeldiagram-bones.hpp"
    #include "modeldiagram-meat.hpp"

    #include "math/distributions.hpp"

    #include "math/lfmcmc.hpp"

    #include "userdata-bones.hpp"
    #include "userdata-meat.hpp"

    #include "seq_processing.hpp"

    #include "database-bones.hpp"
    #include "database-meat.hpp"
    #include "adjlist-bones.hpp"
    #include "adjlist-meat.hpp"

    #include "randgraph.hpp"

    #include "queue-bones.hpp"

    #include "globalevent-bones.hpp"
    #include "globalevent-meat.hpp"

    #include "model-bones.hpp"
    #include "model-meat.hpp"

    #include "viruses-bones.hpp"

    #include "virus-bones.hpp"
    #include "virus-distribute-meat.hpp"
    #include "virus-meat.hpp"
    
    #include "tools-bones.hpp"

    #include "tool-bones.hpp"
    #include "tool-distribute-meat.hpp"
    #include "tool-meat.hpp"

    #include "entity-bones.hpp"
    #include "entity-distribute-meat.hpp"
    #include "entity-meat.hpp"

    #include "entities-bones.hpp"
    
    #include "agent-meat-state.hpp"
    #include "agent-bones.hpp"
    #include "agent-meat.hpp"

    #include "agentssample-bones.hpp"

    #include "models/models.hpp"

}

#endif 

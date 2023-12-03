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

#ifndef EPIWORLD_HPP
#define EPIWORLD_HPP

namespace epiworld {

    #include "config.hpp"
    #include "epiworld-macros.hpp"

    #include "misc.hpp"
    #include "progress.hpp"

    // #include "math/summary-stats.hpp"

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

    #include "globalactions-bones.hpp"
    #include "globalactions-meat.hpp"

    #include "model-bones.hpp"
    #include "model-meat.hpp"

    #include "viruses-bones.hpp"

    #include "virus-bones.hpp"
    #include "virus-meat.hpp"
    
    #include "tools-bones.hpp"

    #include "tool-bones.hpp"
    #include "tool-meat.hpp"

    #include "entity-bones.hpp"
    #include "entity-meat.hpp"

    #include "entities-bones.hpp"
    
    #include "agent-meat-state.hpp"
    #include "agent-bones.hpp"
    #include "agent-meat.hpp"

    #include "agentssample-bones.hpp"

    #include "models/models.hpp"

}

#endif 
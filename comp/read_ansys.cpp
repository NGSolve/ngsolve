#ifdef NGS_PYTHON

#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

#include "../ngstd/python_ngstd.hpp"
#include "meshaccess.hpp"
#include "fespace.hpp"
#include "gridfunction.hpp"
#include <meshing/meshclass.hpp>

namespace ngcomp
{
  using std::string;
  using std::vector;
  using std::map;

  static string Trim (const string & s)
  {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
  }

  static string Upper (string s)
  {
    for (auto & c : s) c = toupper((unsigned char)c);
    return s;
  }

  // collapse internal whitespace to single blanks (for keyword matching)
  static string CollapseWS (const string & s)
  {
    string out;
    bool inws = false;
    for (char c : s)
      {
        if (c == ' ' || c == '\t') { inws = true; continue; }
        if (inws && !out.empty()) out += ' ';
        inws = false;
        out += c;
      }
    return out;
  }

  static vector<string> SplitComma (const string & line)
  {
    vector<string> tok;
    std::stringstream ss(line);
    string item;
    while (std::getline(ss, item, ','))
      tok.push_back(Trim(item));
    return tok;
  }

  static map<string,string> ParseParams (const vector<string> & tok)
  {
    map<string,string> params;
    for (size_t i = 1; i < tok.size(); i++)
      {
        const string & t = tok[i];
        if (t.empty()) continue;
        auto eq = t.find('=');
        if (eq == string::npos)
          params[Upper(Trim(t))] = "";
        else
          params[Upper(Trim(t.substr(0,eq)))] = Trim(t.substr(eq+1));
      }
    return params;
  }

  static string GetParam (const map<string,string> & p, const string & key)
  {
    auto it = p.find(key);
    return (it == p.end()) ? string("") : it->second;
  }


  struct AnsysElement
  {
    netgen::ELEMENT_TYPE type;
    int dim;
    int id;
    vector<int> nodes;        // raw ANSYS node ids, ANSYS node ordering
  };

  struct AnsysMaterial
  {
    double E = 0, nu = 0, rho = 0;
    bool hasE = false, hasNu = false, hasRho = false;
  };

  struct AnsysModel
  {
    vector<int>           node_id;       // ANSYS node ids, in file order
    vector<std::array<double,3>> node_xyz;
    int                   max_node_id = 0;

    vector<AnsysElement>  elements;
    int                   maxdim = 0;

    // named sets (key = display name as first seen)
    map<string, vector<int>> elsets;     // name -> element ids
    map<string, vector<int>> nsets;      // name -> node ids
    // elsets that came from "*ELEMENT, ELSET=" are the bodies/parts (true
    // domains);  standalone "*ELSET" blocks are named selections (for boundaries or other subsets)
    std::set<string> body_elsets;

    // *SOLID/*SHELL SECTION bindings
    vector<std::pair<string,string>> sections;   // (elset name, material name)

    // materials, keyed by display name
    map<string, AnsysMaterial> materials;

    // element id -> index into `elements`
    std::unordered_map<int,size_t> id2elem;

    // *SURFACE, TYPE=ELEMENT : name -> list of (elset-or-elemid token, face id e.g. "S3")
    map<string, vector<std::pair<string,string>>> surfaces;
    // *SURFACE, TYPE=NODE : name -> list of nset-or-nodeid tokens
    map<string, vector<string>> node_surfaces;

    // *CONTACT PAIR / *TIE : (name, interaction, type, slave surface, master surface)
    struct ContactPair { string name, interaction, type, slave, master; };
    vector<ContactPair> contacts;

    // *BOUNDARY : prescribed displacement on a set (dofs 1..3 = x,y,z; 4..6 = rot)
    struct Constraint { string set; vector<int> dofs; double value = 0; };
    vector<Constraint> constraints;

    // *CLOAD : concentrated force,  *DSLOAD : surface pressure,  *DLOAD : body/element load
    struct Load { string kind; string set; int dof = 0; double value = 0; string ltype;
                  double vx = 0, vy = 0, vz = 0; };   // body-force vector (e.g. GRAV)
    vector<Load> loads;
  };


  // ----------------------------------------------------------------------
  //  element type classification
  //
  //  We determine the spatial dimension from the Abaqus type token and the
  //  concrete element type from the actual node count -- this is robust
  //  against the many reduced-integration / hybrid suffixes (C3D8R, C3D10M,
  //  S4R, ...).
  // ----------------------------------------------------------------------
  static int TypeDim (const string & tok)
  {
    if (tok.find("C3D") != string::npos) return 3;   // C3D*, DC3D*, AC3D*
    // 2D continuum / shell / membrane
    static const char* twod[] = { "CPS","CPE","CAX","CPEG","STRI","S","M3D","R3D","SFM" };
    for (auto p : twod)
      if (tok.compare(0, strlen(p), p) == 0) return 2;
    // 1D truss / beam
    if (tok.compare(0,3,"T2D")==0 || tok.compare(0,3,"T3D")==0 ||
        tok.compare(0,1,"B")==0)
      return 1;
    return 0;
  }

  static netgen::ELEMENT_TYPE ClassifyType (const string & rawtok, int nnodes)
  {
    using namespace netgen;
    string tok = Upper(Trim(rawtok));
    int dim = TypeDim(tok);
    switch (dim)
      {
      case 3:
        switch (nnodes)
          {
          case 4:  return TET;
          case 10: return TET10;
          case 8:  return HEX;
          case 20: return HEX20;
          case 6:  return PRISM;
          case 15: return PRISM15;
          case 5:  return PYRAMID;
          }
        break;
      case 2:
        switch (nnodes)
          {
          case 3: return TRIG;
          case 6: return TRIG6;
          case 4: return QUAD;
          case 8: return QUAD8;
          }
        break;
      case 1:
        switch (nnodes)
          {
          case 2: return SEGMENT;
          case 3: return SEGMENT3;
          }
        break;
      }
    throw Exception("ReadAnsysInp: unsupported element type '" + rawtok +
                    "' with " + ToString(nnodes) + " nodes");
  }

  static int ElDim (netgen::ELEMENT_TYPE t)
  {
    using namespace netgen;
    switch (t)
      {
      case SEGMENT: case SEGMENT3: return 1;
      case TRIG: case TRIG6: case QUAD: case QUAD8: return 2;
      default: return 3;
      }
  }

  // permutation: netgen-local node j  <-  abaqus-local node  perm[j].
  static const vector<int> & ReadPerm (netgen::ELEMENT_TYPE t)
  {
    using namespace netgen;
    static const map<netgen::ELEMENT_TYPE, vector<int>> perms = {
      { SEGMENT,  {0,1} },
      { SEGMENT3, {0,1,2} },
      { TRIG,     {0,1,2} },
      { TRIG6,    {0,1,2, 4,5,3} },
      { QUAD,     {0,1,2,3} },
      { QUAD8,    {0,1,2,3, 4,6,7,5} },
      { TET,      {0,1,3,2} },
      { TET10,    {0,1,3,2, 4,7,6,8,5,9} },
      { HEX,      {0,1,2,3,4,5,6,7} },
      { HEX20,    {0,1,2,3,4,5,6,7, 8,10,11,9,12,14,15,13,16,17,18,19} },
      { PRISM,    {1,2,0,4,5,3} },
      { PRISM15,  {1,2,0,4,5,3, 7,6,8,13,14,12,10,9,11} },
      { PYRAMID,  {0,1,2,3,4} },
    };
    auto it = perms.find(t);
    if (it == perms.end())
      throw Exception("ReadAnsysInp: no node ordering for element type " + ToString(int(t)));
    return it->second;
  }

  // number of corner (vertex) nodes of an element type
  static int NumCorners (netgen::ELEMENT_TYPE t)
  {
    using namespace netgen;
    switch (t)
      {
      case SEGMENT: case SEGMENT3:                 return 2;
      case TRIG: case TRIG6:                       return 3;
      case QUAD: case QUAD6: case QUAD8:           return 4;
      case TET: case TET10:                        return 4;
      case PYRAMID: case PYRAMID13:                return 5;
      case PRISM: case PRISM12: case PRISM15:      return 6;
      default:                                     return 8;   // HEX*
      }
  }

  // Abaqus solid face Sk (1-based) -> local corner node indices (0-based, in the
  // element's Abaqus node ordering).
  static vector<int> FaceCornerIndices (netgen::ELEMENT_TYPE t, int sk)
  {
    using namespace netgen;
    switch (t)
      {
      case TET: case TET10:
        switch (sk) { case 1: return {0,1,2}; case 2: return {0,3,1};
                      case 3: return {1,3,2}; case 4: return {2,3,0}; }
        break;
      case HEX: case HEX20:
        switch (sk) { case 1: return {0,1,2,3}; case 2: return {4,5,6,7};
                      case 3: return {0,4,5,1}; case 4: return {1,5,6,2};
                      case 5: return {2,6,7,3}; case 6: return {3,7,4,0}; }
        break;
      case PRISM: case PRISM15:
        switch (sk) { case 1: return {0,1,2}; case 2: return {3,4,5};
                      case 3: return {0,1,4,3}; case 4: return {1,2,5,4}; case 5: return {2,0,3,5}; }
        break;
      case PYRAMID:
        switch (sk) { case 1: return {0,1,2,3}; case 2: return {0,1,4};
                      case 3: return {1,2,4}; case 4: return {2,3,4}; case 5: return {3,0,4}; }
        break;
      default: break;
      }
    return {};
  }


  // ----------------------------------------------------------------------
  //  the Abaqus deck parser
  // ----------------------------------------------------------------------
  static AnsysModel ParseAbaqusInp (const string & filename)
  {
    std::ifstream in(filename);
    if (!in.good())
      throw Exception("ReadAnsysInp: cannot open file '" + filename + "'");

    AnsysModel m;

    enum Mode { M_NONE, M_NODE, M_ELEMENT, M_NSET, M_ELSET, M_ELASTIC, M_DENSITY,
                M_SURFACE, M_SURFACE_NODE, M_CONTACT,
                M_BOUNDARY, M_CLOAD, M_DLOAD, M_DSLOAD, M_SKIP };
    Mode mode = M_NONE;

    // current block context
    netgen::ELEMENT_TYPE cur_eltype = netgen::TET;
    string cur_eltype_token;
    string cur_elem_elset;
    string cur_set_name;
    bool   cur_generate = false;
    string cur_node_nset;       // *NODE, NSET=...
    string cur_material;        // active *MATERIAL
    string cur_surface;         // active *SURFACE name
    AnsysModel::ContactPair cur_contact;   // active *CONTACT PAIR / *TIE template

    auto AddToSet = [](map<string,vector<int>> & sets, const string & name, int id)
    {
      if (!name.empty()) sets[name].push_back(id);
    };

    string pending;
    bool has_pending = false;
    auto NextLine = [&](string & out) -> bool
    {
      if (has_pending) { out = pending; has_pending = false; return true; }
      return (bool)std::getline(in, out);
    };

    string raw;
    while (NextLine(raw))
      {
        string line = Trim(raw);
        if (line.empty()) continue;
        if (line.size() >= 2 && line[0] == '*' && line[1] == '*') continue; // comment

        if (line[0] == '*')
          {
            // ---- keyword line ----
            auto tok = SplitComma(line);
            string kw = Upper(CollapseWS(tok[0]));
            auto params = ParseParams(tok);

            if (kw == "*NODE")
              {
                mode = M_NODE;
                cur_node_nset = GetParam(params, "NSET");
              }
            else if (kw == "*ELEMENT")
              {
                mode = M_ELEMENT;
                cur_eltype_token = GetParam(params, "TYPE");
                cur_elem_elset   = GetParam(params, "ELSET");
                if (!cur_elem_elset.empty()) m.body_elsets.insert(cur_elem_elset);
                // concrete type resolved per data-line from node count
              }
            else if (kw == "*NSET")
              {
                mode = M_NSET;
                cur_set_name = GetParam(params, "NSET");
                cur_generate = params.count("GENERATE") > 0;
                if (!cur_set_name.empty()) m.nsets[cur_set_name];   // ensure exists
              }
            else if (kw == "*ELSET")
              {
                mode = M_ELSET;
                cur_set_name = GetParam(params, "ELSET");
                cur_generate = params.count("GENERATE") > 0;
                if (!cur_set_name.empty()) m.elsets[cur_set_name];
              }
            else if (kw == "*SOLID SECTION" || kw == "*SHELL SECTION")
              {
                string es = GetParam(params, "ELSET");
                string mt = GetParam(params, "MATERIAL");
                if (!es.empty()) m.sections.push_back({es, mt});
                mode = M_SKIP;     // ignore the (thickness) data line
              }
            else if (kw == "*MATERIAL")
              {
                cur_material = GetParam(params, "NAME");
                if (!cur_material.empty()) m.materials[cur_material];
                mode = M_NONE;
              }
            else if (kw == "*ELASTIC")
              mode = M_ELASTIC;
            else if (kw == "*DENSITY")
              mode = M_DENSITY;
            else if (kw == "*SURFACE")
              {
                cur_surface = GetParam(params, "NAME");
                string stype = Upper(GetParam(params, "TYPE"));
                if (stype == "NODE")
                  { mode = M_SURFACE_NODE; if (!cur_surface.empty()) m.node_surfaces[cur_surface]; }
                else
                  { mode = M_SURFACE;      if (!cur_surface.empty()) m.surfaces[cur_surface]; }
              }
            else if (kw == "*CONTACT PAIR" || kw == "*TIE")
              {
                cur_contact = AnsysModel::ContactPair{};
                cur_contact.name        = GetParam(params, "NAME");
                cur_contact.interaction = GetParam(params, "INTERACTION");
                cur_contact.type        = (kw == "*TIE") ? string("TIE") : GetParam(params, "TYPE");
                mode = M_CONTACT;
              }
            else if (kw == "*BOUNDARY")  mode = M_BOUNDARY;
            else if (kw == "*CLOAD")     mode = M_CLOAD;
            else if (kw == "*DLOAD")     mode = M_DLOAD;
            else if (kw == "*DSLOAD")    mode = M_DSLOAD;
            else
              mode = M_SKIP;       // *HEADING, *PART, *STEP, *SURFACE INTERACTION, ...
            continue;
          }

        // ---- data line; handle Abaqus continuation (trailing comma) ----
        while (!line.empty() && line.back() == ',')
          {
            string next;
            if (!NextLine(next)) break;
            string tn = Trim(next);
            if (!tn.empty() && tn[0] == '*')      // keyword/comment -> not a continuation
              { pending = next; has_pending = true; break; }
            line += tn;
          }
        auto tok = SplitComma(line);
        // drop trailing empty token from a trailing comma
        while (!tok.empty() && tok.back().empty()) tok.pop_back();
        if (tok.empty()) continue;

        switch (mode)
          {
          case M_NODE:
            {
              int id = std::stoi(tok[0]);
              std::array<double,3> p = {0,0,0};
              for (int k = 0; k < 3 && k+1 < (int)tok.size(); k++)
                p[k] = std::stod(tok[k+1]);
              m.node_id.push_back(id);
              m.node_xyz.push_back(p);
              m.max_node_id = max(m.max_node_id, id);
              AddToSet(m.nsets, cur_node_nset, id);
              break;
            }
          case M_ELEMENT:
            {
              AnsysElement el;
              el.id = std::stoi(tok[0]);
              for (size_t k = 1; k < tok.size(); k++)
                if (!tok[k].empty()) el.nodes.push_back(std::stoi(tok[k]));
              el.type = ClassifyType(cur_eltype_token, el.nodes.size());
              el.dim  = ElDim(el.type);
              m.maxdim = max(m.maxdim, el.dim);
              m.id2elem[el.id] = m.elements.size();
              m.elements.push_back(std::move(el));
              AddToSet(m.elsets, cur_elem_elset, m.elements.back().id);
              break;
            }
          case M_NSET:
          case M_ELSET:
            {
              auto & sets = (mode == M_NSET) ? m.nsets : m.elsets;
              auto & vec  = sets[cur_set_name];
              if (cur_generate && tok.size() >= 2)
                {
                  int a = std::stoi(tok[0]), b = std::stoi(tok[1]);
                  int s = (tok.size() >= 3 && !tok[2].empty()) ? std::stoi(tok[2]) : 1;
                  if (s == 0) s = 1;
                  for (int v = a; v <= b; v += s) vec.push_back(v);
                }
              else
                for (auto & t : tok)
                  {
                    if (t.empty()) continue;
                    try { vec.push_back(std::stoi(t)); }
                    catch (...) { /* set-of-sets reference: phase 2 */ }
                  }
              break;
            }
          case M_ELASTIC:
            {
              auto & mat = m.materials[cur_material];
              if (tok.size() >= 1 && !tok[0].empty()) { mat.E  = std::stod(tok[0]); mat.hasE  = true; }
              if (tok.size() >= 2 && !tok[1].empty()) { mat.nu = std::stod(tok[1]); mat.hasNu = true; }
              mode = M_SKIP;   // ignore further temperature lines (phase 1)
              break;
            }
          case M_DENSITY:
            {
              auto & mat = m.materials[cur_material];
              if (tok.size() >= 1 && !tok[0].empty()) { mat.rho = std::stod(tok[0]); mat.hasRho = true; }
              mode = M_SKIP;
              break;
            }
          case M_SURFACE:        // data: elset-or-elemid, face-id (e.g. S3)
            {
              if (!cur_surface.empty() && tok.size() >= 1)
                {
                  string face = (tok.size() >= 2) ? Upper(tok[1]) : "";
                  m.surfaces[cur_surface].push_back({ tok[0], face });
                }
              break;
            }
          case M_SURFACE_NODE:   // data: nset-or-nodeid
            {
              if (!cur_surface.empty() && tok.size() >= 1 && !tok[0].empty())
                m.node_surfaces[cur_surface].push_back(tok[0]);
              break;
            }
          case M_CONTACT:        // data: slave-surface, master-surface
            {
              if (tok.size() >= 2 && !tok[0].empty() && !tok[1].empty())
                {
                  auto c = cur_contact;
                  c.slave = tok[0];
                  c.master = tok[1];
                  m.contacts.push_back(c);
                }
              break;
            }
          case M_BOUNDARY:       // data: set, dof1[, dof2[, value]]  or  set, TYPE
            {
              if (tok.size() >= 2 && !tok[0].empty())
                {
                  AnsysModel::Constraint cn;
                  cn.set = tok[0];
                  string t1 = Upper(tok[1]);
                  // named conditions
                  if (t1 == "ENCASTRE")      cn.dofs = {1,2,3,4,5,6};
                  else if (t1 == "PINNED")   cn.dofs = {1,2,3};
                  else if (t1 == "XSYMM")    cn.dofs = {1,5,6};
                  else if (t1 == "YSYMM")    cn.dofs = {2,4,6};
                  else if (t1 == "ZSYMM")    cn.dofs = {3,4,5};
                  else
                    {
                      int d1 = std::stoi(tok[1]);
                      int d2 = (tok.size() >= 3 && !tok[2].empty()) ? std::stoi(tok[2]) : d1;
                      for (int d = d1; d <= d2; d++) cn.dofs.push_back(d);
                      if (tok.size() >= 4 && !tok[3].empty()) cn.value = std::stod(tok[3]);
                    }
                  m.constraints.push_back(std::move(cn));
                }
              break;
            }
          case M_CLOAD:          // data: set, dof, magnitude
            {
              if (tok.size() >= 3 && !tok[0].empty())
                m.loads.push_back({ "force", tok[0], std::stoi(tok[1]), std::stod(tok[2]), "" });
              break;
            }
          case M_DSLOAD:         // data: surface, Ptype, magnitude
            {
              if (tok.size() >= 3 && !tok[0].empty())
                m.loads.push_back({ "pressure", tok[0], 0, std::stod(tok[2]), Upper(tok[1]) });
              break;
            }
          case M_DLOAD:          // data: elset, Ptype, magnitude[, dir...]
            {
              if (tok.size() >= 2 && !tok[0].empty())
                {
                  string lt = Upper(tok[1]);
                  if (lt.size() >= 2 && lt[0] == 'P' && isdigit((unsigned char)lt[1]))
                    {
                      // element-face pressure Pn -> synthesise a *SURFACE so it
                      // resolves to a boundary region like *DSLOAD
                      string synth = tok[0] + "_" + lt;
                      m.surfaces[synth].push_back({ tok[0], "S" + lt.substr(1) });
                      double v = (tok.size() >= 3 && !tok[2].empty()) ? std::stod(tok[2]) : 0;
                      m.loads.push_back({ "pressure", synth, 0, v, lt });
                    }
                  else
                    {
                      AnsysModel::Load l; l.kind = "body"; l.set = tok[0]; l.ltype = lt;
                      if (lt == "GRAV" && tok.size() >= 6)
                        { double g = std::stod(tok[2]); l.value = g;
                          l.vx = g*std::stod(tok[3]); l.vy = g*std::stod(tok[4]); l.vz = g*std::stod(tok[5]); }
                      else if (tok.size() >= 3 && !tok[2].empty())
                        { try { l.value = std::stod(tok[2]); } catch (...) {} }
                      m.loads.push_back(std::move(l));
                    }
                }
              break;
            }
          default:
            break;   // M_NONE / M_SKIP
          }
      }

    return m;
  }


  // ----------------------------------------------------------------------
  //  build a netgen mesh from the parsed model
  //
  //  Out parameters carry the bookkeeping needed to build Region objects:
  //    vol_regions[name]  : 0-based VOL region numbers carrying that elset
  //    bnd_regions[name]  : 0-based BND region numbers carrying that elset
  //    mat_cells[matname] : 0-based VOL region numbers using that material
  //    id2vertex[ansysid] : 0-based ngsolve vertex number (-1 if absent)
  // ----------------------------------------------------------------------
  static shared_ptr<netgen::Mesh>
  BuildAnsysMesh (const AnsysModel & m, double scale,
                  map<string,vector<int>> & vol_regions,
                  map<string,vector<int>> & bnd_regions,
                  map<string,vector<int>> & bbnd_regions,
                  map<string,vector<int>> & mat_cells,
                  map<string,vector<int>> & element_sets,
                  vector<int> & id2vertex)
  {
    using namespace netgen;

    const int meshdim = m.maxdim;
    if (meshdim != 2 && meshdim != 3)
      throw Exception("ReadAnsysInp: need a 2D or 3D mesh "
                      "(found max element dimension " + ToString(meshdim) + ")");

    auto mesh = make_shared<Mesh>();
    mesh->SetDimension(meshdim);

    // ---- points + ANSYS-id -> netgen point remap ----
    vector<int> id2pi(m.max_node_id+1, -1);
    id2vertex.assign(m.max_node_id+1, -1);
    for (size_t i = 0; i < m.node_id.size(); i++)
      {
        const auto & c = m.node_xyz[i];
        PointIndex pi = mesh->AddPoint(Point3d(scale*c[0], scale*c[1], scale*c[2]));
        id2pi[m.node_id[i]]     = (int)pi;
        id2vertex[m.node_id[i]] = (int)pi - PointIndex::BASE;   // 0-based vertex nr
      }

    // ---- reverse membership: element id -> set names containing it ----
    std::unordered_map<int, vector<string>> elem2sets;     // all elsets
    std::unordered_map<int, vector<string>> elem2bodies;   // only *ELEMENT (body) elsets
    for (auto & [name, ids] : m.elsets)
      {
        bool body = m.body_elsets.count(name) > 0;
        for (int id : ids)
          {
            elem2sets[id].push_back(name);
            if (body) elem2bodies[id].push_back(name);
          }
      }

    // element id -> material (from any *SECTION, resolving its elset)
    std::unordered_map<int, string> elem_material;
    for (auto & [es, mat] : m.sections)
      {
        if (mat.empty()) continue;
        auto it = m.elsets.find(es);
        if (it == m.elsets.end()) continue;
        for (int id : it->second) elem_material.emplace(id, mat);
      }
    auto MaterialOf = [&](int elemid) -> string
    {
      auto it = elem_material.find(elemid);
      return (it != elem_material.end()) ? it->second : "";
    };

    // full label (all sets) and body label (only *ELEMENT/part elsets)
    auto SortedLabel = [&](int elemid) -> vector<string>
    {
      vector<string> lab;
      auto it = elem2sets.find(elemid);
      if (it != elem2sets.end()) lab = it->second;
      std::sort(lab.begin(), lab.end());
      return lab;
    };
    auto BodyLabel = [&](int elemid) -> vector<string>
    {
      vector<string> lab;
      auto it = elem2bodies.find(elemid);
      if (it != elem2bodies.end()) lab = it->second;
      std::sort(lab.begin(), lab.end());
      return lab;
    };

    // "size" of a named group: how many elements / faces it contains.  Bigger =
    // more general (e.g. a whole body), smaller = more specific (e.g. a face set).
    auto GroupSize = [&](const string & name) -> size_t
    {
      auto e = m.elsets.find(name);   if (e != m.elsets.end())   return e->second.size();
      auto s = m.surfaces.find(name); if (s != m.surfaces.end()) return s->second.size();
      return 0;
    };

    // compose a region name from its set-labels, ordered most-general first
    // (largest group) to most-specific last, so e.g. "Solid_part-2+contact_S1".
    auto JoinLabel = [&](vector<string> label, const string & dflt) -> string
    {
      if (label.empty()) return dflt;
      std::sort(label.begin(), label.end(), [&](const string & a, const string & b)
                { size_t sa = GroupSize(a), sb = GroupSize(b);
                  return (sa != sb) ? (sa > sb) : (a < b); });
      string name = label[0];
      for (size_t k = 1; k < label.size(); k++) name += "+" + label[k];
      return name;
    };

    auto SetPNums = [&](auto & el, const AnsysElement & ae)
    {
      const auto & perm = ReadPerm(ae.type);
      for (size_t j = 0; j < perm.size(); j++)
        {
          int nid = ae.nodes[perm[j]];
          int pid = (nid >= 1 && nid < (int)id2pi.size()) ? id2pi[nid] : -1;
          if (pid < 0) throw Exception("ReadAnsysInp: element " + ToString(ae.id) +
                                       " references undefined node " + ToString(nid));
          el[j] = PointIndex(pid);
        }
    };

    // ---- generic finest-partition over the named sets containing each element
    //      of a given spatial dimension; calls onNewRegion(regnr,name) the first
    //      time a partition cell is seen and onElement(ae,regnr) for every element.
    auto Partition = [&](int eldim, bool isMaterial,
                         const std::function<vector<string>(int)> & labelOf,
                         map<string,vector<int>> & regions_out,
                         const std::function<void(int,const string&)> & onNewRegion,
                         const std::function<void(const AnsysElement&,int)> & onElement)
    {
      map<vector<string>, int> label2reg;
      for (auto & ae : m.elements)
        {
          if (ae.dim != eldim) continue;
          auto label = labelOf(ae.id);
          int regnr;
          auto it = label2reg.find(label);
          if (it == label2reg.end())
            {
              regnr = label2reg.size();
              label2reg[label] = regnr;
              string mat  = isMaterial ? MaterialOf(ae.id) : "";
              string name = mat.empty() ? JoinLabel(label, "default") : mat;
              for (auto & l : label) regions_out[l].push_back(regnr);
              if (isMaterial && !mat.empty()) mat_cells[mat].push_back(regnr);
              onNewRegion(regnr, name);
            }
          else
            regnr = it->second;
          onElement(ae, regnr);
        }
      return (int)label2reg.size();
    };

    // ====================  domain (codim 0) elements  ====================
    // partition by body (part) only -> clean "true" domains, named by material
    std::unordered_map<int,int> elemid2elnr;   // ANSYS elem id -> ngsolve VOL element nr
    int nvol_added = 0;
    Partition(meshdim, /*isMaterial*/true, BodyLabel, vol_regions,
      [&](int regnr, const string & name)
      {
        mesh->SetMaterial(regnr+1, name);                       // 1-based
        if (meshdim == 2)
          {
            int fdnr = mesh->AddFaceDescriptor(FaceDescriptor(regnr+1, regnr+1, 0, 0));
            mesh->GetFaceDescriptor(fdnr).SetBCProperty(regnr+1);
          }
      },
      [&](const AnsysElement & ae, int regnr)
      {
        elemid2elnr[ae.id] = nvol_added++;                 // 0-based ngsolve VOL element nr
        if (meshdim == 3)
          { Element el(ae.type);   SetPNums(el, ae); el.SetIndex(regnr+1); mesh->AddVolumeElement(el); }
        else
          { Element2d el(ae.type); SetPNums(el, ae); el.SetIndex(regnr+1); mesh->AddSurfaceElement(el); }
      });

    // standalone *ELSET selections (not bodies) -> element-number lists, so the
    // many auto-generated named selections do not pollute the domain partition
    for (auto & [name, ids] : m.elsets)
      if (!m.body_elsets.count(name))
        {
          vector<int> els;
          for (int id : ids)
            { auto it = elemid2elnr.find(id); if (it != elemid2elnr.end()) els.push_back(it->second); }
          if (!els.empty()) element_sets[name] = std::move(els);
        }

    // vertices must be counted before FindOpenElements/FindOpenSegments can
    // identify the boundary correctly
    mesh->ComputeNVertices();

    bool have_surface = false;

    if (meshdim == 3)
      {
        // ===========  edge (codim 2) = 1D elements (line loads, edge BCs) ======
        Partition(1, /*isMaterial*/false, BodyLabel, bbnd_regions,
          [&](int regnr, const string & name) { mesh->SetCD2Name(regnr+1, name); },  // 1-based
          [&](const AnsysElement & ae, int regnr)
          { Segment seg; SetPNums(seg, ae); seg.SetIndex(regnr+1); mesh->AddSegment(seg); });

        // ===========  boundary (codim 1) =====================================
        // The actual mesh boundary is the set of volume faces occurring exactly
        // once.  We enumerate them with Element::GetFace (correct TRIG6/QUAD8 for
        // high-order elements), then *name* each one by matching its corner
        // vertices against named sources -- *SURFACE element faces and explicit
        // shell elements -- so contact / load surfaces become boundary regions.
        // (FindOpenElements is avoided: it mis-builds high-order faces as quads.)
        auto AnsysKey = [&](const vector<int> & nids) -> vector<int>
        {
          vector<int> key;
          for (int nid : nids)
            {
              int pid = (nid >= 1 && nid < (int)id2pi.size()) ? id2pi[nid] : -1;
              if (pid < 0) return {};
              key.push_back(pid);
            }
          std::sort(key.begin(), key.end());
          return key;
        };
        auto CornerKey = [&](const Element2d & f) -> vector<int>
        {
          int nc = (f.GetType()==QUAD || f.GetType()==QUAD6 || f.GetType()==QUAD8) ? 4 : 3;
          vector<int> key;
          for (int c = 0; c < nc; c++) key.push_back((int)f[c]);
          std::sort(key.begin(), key.end());
          return key;
        };

        // named facets: corner key -> set of region names
        std::map<vector<int>, std::set<string>> namedlabel;
        // (a) *SURFACE element faces
        for (auto & [sname, entries] : m.surfaces)
          for (auto & [token, faceid] : entries)
            {
              vector<int> elemids;
              auto se = m.elsets.find(token);
              if (se != m.elsets.end()) elemids = se->second;
              else { try { elemids.push_back(std::stoi(token)); } catch (...) {} }
              bool allcorners = !(faceid.size() >= 2 && faceid[0] == 'S'
                                  && isdigit((unsigned char)faceid[1]));
              int sk = allcorners ? 0 : std::stoi(faceid.substr(1));
              for (int eid : elemids)
                {
                  auto it = m.id2elem.find(eid);
                  if (it == m.id2elem.end()) continue;
                  const AnsysElement & ae = m.elements[it->second];
                  vector<int> loc = allcorners
                    ? [&]{ vector<int> v; for (int k=0;k<NumCorners(ae.type);k++) v.push_back(k); return v; }()
                    : FaceCornerIndices(ae.type, sk);
                  vector<int> corners;
                  for (int l : loc) if (l < (int)ae.nodes.size()) corners.push_back(ae.nodes[l]);
                  auto key = AnsysKey(corners);
                  if (!key.empty()) namedlabel[key].insert(sname);
                }
            }
        // (b) explicit shell / 2D elements
        for (auto & ae : m.elements)
          if (ae.dim == 2)
            {
              int nc = NumCorners(ae.type);
              vector<int> corners(ae.nodes.begin(),
                                  ae.nodes.begin() + min((size_t)nc, ae.nodes.size()));
              auto key = AnsysKey(corners);
              if (key.empty()) continue;
              for (auto & l : SortedLabel(ae.id)) namedlabel[key].insert(l);
            }

        // enumerate volume faces, counting by corner key
        struct FaceRec { Element2d face; ElementIndex ei; int cnt; };
        std::map<vector<int>, FaceRec> faces;
        for (ElementIndex ei : mesh->VolumeElements().Range())
          {
            const Element & el = (*mesh)[ei];
            for (int j = 1; j <= el.GetNFaces(); j++)
              {
                Element2d f;
                el.GetFace(j, f);
                auto key = CornerKey(f);
                auto it = faces.find(key);
                if (it == faces.end()) faces[key] = FaceRec{ f, ei, 1 };
                else                   it->second.cnt++;
              }
          }

        // orient a boundary face outward (normal away from element centroid --
        // robust to GetFace's winding, which differs for high-order elements)
        auto OrientOutward = [&](Element2d & f, ElementIndex ei)
        {
          const Element & el = (*mesh)[ei];
          double ex=0, ey=0, ez=0; int np = el.GetNP();
          for (int k = 0; k < np; k++)
            { auto p = mesh->Point(el[k]); ex+=p(0); ey+=p(1); ez+=p(2); }
          ex/=np; ey/=np; ez/=np;
          auto a = mesh->Point(f[0]), b = mesh->Point(f[1]), c = mesh->Point(f[2]);
          double ux=b(0)-a(0), uy=b(1)-a(1), uz=b(2)-a(2);
          double vx=c(0)-a(0), vy=c(1)-a(1), vz=c(2)-a(2);
          double nx=uy*vz-uz*vy, ny=uz*vx-ux*vz, nz=ux*vy-uy*vx;
          double fx=(a(0)+b(0)+c(0))/3-ex, fy=(a(1)+b(1)+c(1))/3-ey, fz=(a(2)+b(2)+c(2))/3-ez;
          if (nx*fx+ny*fy+nz*fz < 0) f.Invert();
        };

        // boundary region per (sorted) label;  one FaceDescriptor per (region,domain)
        std::map<vector<string>, int> label2reg;
        std::map<std::pair<int,int>, int> regdom2fd;
        auto GetRegion = [&](const vector<string> & label) -> int
        {
          auto it = label2reg.find(label);
          if (it != label2reg.end()) return it->second;
          int regnr = label2reg.size();
          label2reg[label] = regnr;
          mesh->SetBCName(regnr, JoinLabel(label, "default"));
          if (label.empty()) bnd_regions["default"].push_back(regnr);
          for (auto & l : label) bnd_regions[l].push_back(regnr);
          return regnr;
        };
        auto AddBndFace = [&](Element2d & f, int regnr, int dom)
        {
          auto k = std::make_pair(regnr, dom);
          int fdnr;
          auto it = regdom2fd.find(k);
          if (it == regdom2fd.end())
            {
              fdnr = mesh->AddFaceDescriptor(FaceDescriptor(regnr+1, dom, 0, 0));
              mesh->GetFaceDescriptor(fdnr).SetBCProperty(regnr+1);
              regdom2fd[k] = fdnr;
            }
          else fdnr = it->second;
          f.SetIndex(fdnr);
          mesh->AddSurfaceElement(f);
          have_surface = true;
        };
        // node-based references (node sets / node *SURFACEs) that boundary
        // conditions point at: a boundary facet belongs to such a set iff *all*
        // its corner nodes are in the set.  Only sets actually referenced by a
        // BC/load/node-surface are resolved, to keep the boundary partition small.
        std::unordered_map<int, vector<string>> pi2nodenames;   // netgen pi -> covering set names
        {
          std::set<string> refs;
          for (auto & c : m.constraints)     refs.insert(c.set);
          for (auto & l : m.loads)           refs.insert(l.set);
          for (auto & [n, t] : m.node_surfaces) refs.insert(n);
          auto addNodes = [&](const string & name, const vector<int> & ids)
          {
            for (int id : ids)
              { int pid = (id>=1 && id<(int)id2pi.size()) ? id2pi[id] : -1;
                if (pid >= 0) pi2nodenames[pid].push_back(name); }
          };
          for (auto & name : refs)
            {
              auto nit = m.nsets.find(name);
              if (nit != m.nsets.end()) { addNodes(name, nit->second); continue; }
              auto sit = m.node_surfaces.find(name);
              if (sit != m.node_surfaces.end())
                {
                  vector<int> ids;
                  for (auto & t : sit->second)
                    { auto e = m.nsets.find(t);
                      if (e != m.nsets.end()) ids.insert(ids.end(), e->second.begin(), e->second.end());
                      else { try { ids.push_back(std::stoi(t)); } catch (...) {} } }
                  addNodes(name, ids);
                }
            }
        }
        // full (sorted) label of a boundary facet: *SURFACE/shell names + covering node sets
        auto FacetLabel = [&](const vector<int> & key, const Element2d & f) -> vector<string>
        {
          std::set<string> lab;
          auto nl = namedlabel.find(key);
          if (nl != namedlabel.end()) lab.insert(nl->second.begin(), nl->second.end());
          if (!pi2nodenames.empty())
            {
              int nc = (f.GetType()==QUAD || f.GetType()==QUAD6 || f.GetType()==QUAD8) ? 4 : 3;
              std::map<string,int> cnt;
              for (int c = 0; c < nc; c++)
                { auto it = pi2nodenames.find((int)f[c]);
                  if (it != pi2nodenames.end()) for (auto & nm : it->second) cnt[nm]++; }
              for (auto & [nm, k] : cnt) if (k == nc) lab.insert(nm);   // all corners covered
            }
          return vector<string>(lab.begin(), lab.end());
        };

        std::set<vector<int>> usedkeys;
        for (auto & [key, fr] : faces)
          {
            if (fr.cnt != 1) continue;                    // interior face
            usedkeys.insert(key);
            int regnr = GetRegion(FacetLabel(key, fr.face));
            Element2d f = fr.face;
            OrientOutward(f, fr.ei);
            AddBndFace(f, regnr, (*mesh)[fr.ei].GetIndex());
          }
        // standalone shells (no adjacent volume face -- e.g. shell-only models)
        for (auto & ae : m.elements)
          if (ae.dim == 2)
            {
              int nc = NumCorners(ae.type);
              vector<int> corners(ae.nodes.begin(),
                                  ae.nodes.begin() + min((size_t)nc, ae.nodes.size()));
              auto key = AnsysKey(corners);
              if (key.empty() || usedkeys.count(key)) continue;
              usedkeys.insert(key);
              Element2d el(ae.type); SetPNums(el, ae);
              int regnr = GetRegion(FacetLabel(key, el));
              AddBndFace(el, regnr, 0);
            }
      }
    else // meshdim == 2: boundary (codim 1) = segments
      {
        std::map<std::pair<int,int>, vector<string>> namededge;   // sorted ansys node pair -> label
        for (auto & ae : m.elements)
          if (ae.dim == 1 && ae.nodes.size() >= 2)
            {
              int a = ae.nodes[0], b = ae.nodes[1];
              namededge[{min(a,b), max(a,b)}] = SortedLabel(ae.id);
            }

        // reverse map netgen point -> ansys node id
        vector<int> pi2id(mesh->GetNP()+2, -1);
        for (size_t i = 0; i < m.node_id.size(); i++)
          pi2id[id2pi[m.node_id[i]]] = m.node_id[i];

        mesh->CalcSurfacesOfNode();          // prerequisite for FindOpenSegments
        mesh->FindOpenSegments();
        int nopen = mesh->GetNOpenSegments();
        map<vector<string>, int> label2reg;
        for (int i = 1; i <= nopen; i++)
          {
            Segment seg = mesh->GetOpenSegment(i);
            int a = pi2id[(int)seg[0]], b = pi2id[(int)seg[1]];
            vector<string> label;
            if (a > 0 && b > 0)
              {
                auto it = namededge.find({min(a,b), max(a,b)});
                if (it != namededge.end()) label = it->second;
              }
            int regnr;
            auto it = label2reg.find(label);
            if (it == label2reg.end())
              {
                regnr = label2reg.size();
                label2reg[label] = regnr;
                string name = JoinLabel(label, "default");
                mesh->SetBCName(regnr, name);                   // 0-based
                EdgeDescriptor ed; ed.SetName(name);
                mesh->AddEdgeDescriptor(ed);                    // -> 1-based index regnr+1
                for (auto & l : label) bnd_regions[l].push_back(regnr);
                if (label.empty()) bnd_regions["default"].push_back(regnr);
              }
            else regnr = it->second;
            seg.SetIndex(regnr+1);
            mesh->AddSegment(seg);
            have_surface = true;
          }
      }

    // ---- finalize ----
    if (have_surface) mesh->RebuildSurfaceElementLists();
    Point3d pmin, pmax;
    mesh->GetBox(pmin, pmax);
    mesh->UpdateTopology();

    return mesh;
  }


} // namespace ngcomp


using namespace ngcomp;

void ExportReadAnsys (py::module & m)
{
    m.def("ReadAnsysInp", [](const string & filename, double scale) -> py::dict
    {
      AnsysModel model = ParseAbaqusInp(filename);

      map<string,vector<int>> vol_regions, bnd_regions, bbnd_regions, mat_cells, element_sets;
      vector<int> id2vertex;
      auto ngmesh = BuildAnsysMesh(model, scale, vol_regions, bnd_regions,
                                   bbnd_regions, mat_cells, element_sets, id2vertex);

      auto ma = make_shared<MeshAccess>(ngmesh);
      ngmesh->updateSignal.Connect(ma.get(), [p=ma.get()] { p->UpdateBuffers(); });

      const int nvol  = ma->GetNRegions(VOL);
      const int nbnd  = ma->GetNRegions(BND);
      const int nbbnd = ma->GetNRegions(BBND);

      auto MakeRegion = [&](VorB vb, const vector<int> & nrs, int nregions)
      {
        BitArray ba(nregions);
        ba.Clear();
        for (int r : nrs) if (r >= 0 && r < nregions) ba.SetBit(r);
        return Region(ma, vb, ba);
      };

      py::dict res;
      res["mesh"] = py::cast(ma);

      // material property data
      py::dict mats;
      for (auto & [name, mat] : model.materials)
        {
          py::dict d;
          if (mat.hasE)   d["E"]   = mat.E;
          if (mat.hasNu)  d["nu"]  = mat.nu;
          if (mat.hasRho) d["rho"] = mat.rho;
          mats[py::str(name)] = d;
        }
      res["materials"] = mats;

      // "true" domains: one Region(VOL) per body/part (overlap-safe)
      py::dict domains;
      for (auto & [name, nrs] : vol_regions)
        domains[py::str(name)] = py::cast(MakeRegion(VOL, nrs, nvol));
      res["domains"] = domains;

      py::dict material_regions;
      for (auto & [name, nrs] : mat_cells)
        material_regions[py::str(name)] = py::cast(MakeRegion(VOL, nrs, nvol));
      res["material_regions"] = material_regions;

      // standalone *ELSET named selections (not bodies): element-number lists,
      // kept out of the domain partition so it stays clean
      py::dict elsets;
      for (auto & [name, els] : element_sets)
        {
          py::list l;
          for (int e : els) l.append(e);
          elsets[py::str(name)] = l;
        }
      res["element_sets"] = elsets;

      py::dict boundaries;
      for (auto & [name, nrs] : bnd_regions)
        boundaries[py::str(name)] = py::cast(MakeRegion(BND, nrs, nbnd));
      res["boundaries"] = boundaries;

      // contact pairs / ties: each references two *SURFACE names, which are
      // boundary regions in `boundaries`
      py::list contacts;
      for (auto & c : model.contacts)
        {
          auto reg = [&](const string & sname) -> py::object
          {
            auto it = bnd_regions.find(sname);
            if (it == bnd_regions.end()) return py::none();
            return py::cast(MakeRegion(BND, it->second, nbnd));
          };
          py::dict d;
          d["name"]        = c.name;
          d["interaction"] = c.interaction;
          d["type"]        = c.type;
          d["slave"]       = reg(c.slave);
          d["master"]      = reg(c.master);
          d["slave_name"]  = c.slave;
          d["master_name"] = c.master;
          contacts.append(d);
        }
      res["contact_pairs"] = contacts;

      // codim-2 edge regions (3D only; e.g. for line loads)
      py::dict edges;
      for (auto & [name, nrs] : bbnd_regions)
        edges[py::str(name)] = py::cast(MakeRegion(BBND, nrs, nbbnd));
      res["edges"] = edges;

      // node sets -> 0-based ngsolve vertex numbers
      py::dict nodesets;
      for (auto & [name, ids] : model.nsets)
        {
          py::list verts;
          for (int id : ids)
            if (id >= 0 && id < (int)id2vertex.size() && id2vertex[id] >= 0)
              verts.append(id2vertex[id]);
          nodesets[py::str(name)] = verts;
        }
      res["node_sets"] = nodesets;

      // resolve a set name to the most specific region we have for it
      auto Resolve = [&](const string & name, py::dict & d)
      {
        d["name"] = name;
        auto b = bnd_regions.find(name);
        if (b != bnd_regions.end())
          { d["region"] = py::cast(MakeRegion(BND, b->second, nbnd)); d["vb"] = "BND"; return; }
        auto v = vol_regions.find(name);
        if (v != vol_regions.end())
          { d["region"] = py::cast(MakeRegion(VOL, v->second, nvol)); d["vb"] = "VOL"; return; }
        auto e = bbnd_regions.find(name);
        if (e != bbnd_regions.end())
          { d["region"] = py::cast(MakeRegion(BBND, e->second, nbbnd)); d["vb"] = "BBND"; return; }
        auto n = model.nsets.find(name);
        if (n != model.nsets.end())
          {
            py::list verts;
            for (int id : n->second)
              if (id >= 0 && id < (int)id2vertex.size() && id2vertex[id] >= 0)
                verts.append(id2vertex[id]);
            d["region"] = py::none(); d["vb"] = "nodes"; d["nodes"] = verts; return;
          }
        d["region"] = py::none(); d["vb"] = py::none();
      };

      // *BOUNDARY -> prescribed-displacement constraints (fixed supports etc.)
      py::list constraints;
      for (auto & c : model.constraints)
        {
          py::dict d;
          Resolve(c.set, d);
          py::list dofs;
          for (int dd : c.dofs) dofs.append(dd);
          d["dofs"]  = dofs;          // 1,2,3 = x,y,z ; 4,5,6 = rotations
          d["value"] = c.value;
          constraints.append(d);
        }
      res["constraints"] = constraints;

      // *CLOAD / *DSLOAD / *DLOAD -> loads
      auto Const = [](double v) -> shared_ptr<CoefficientFunction>
      { return make_shared<ConstantCoefficientFunction>(v); };   // real-valued
      auto VecCF = [&](double a, double b, double c)
      {
        Array<shared_ptr<CoefficientFunction>> comps;
        comps.Append(Const(a)); comps.Append(Const(b)); comps.Append(Const(c));
        return MakeVectorialCoefficientFunction(std::move(comps));
      };

      py::list loads;
      Array<pair<variant<string,Region>, shared_ptr<CoefficientFunction>>> pvals, bvals;
      for (auto & l : model.loads)
        {
          py::dict d;
          Resolve(l.set, d);
          d["kind"]  = l.kind;        // "force" | "pressure" | "body"
          d["value"] = l.value;
          if (l.kind == "force") d["dof"] = l.dof;       // 1,2,3 = x,y,z
          if (!l.ltype.empty())  d["ltype"] = l.ltype;   // e.g. "P", "GRAV"
          if (l.kind == "body")
            {
              d["vector"] = py::make_tuple(l.vx, l.vy, l.vz);
              auto v = vol_regions.find(l.set);
              if (v != vol_regions.end())
                bvals.Append({ MakeRegion(VOL, v->second, nvol), VecCF(l.vx, l.vy, l.vz) });
            }
          else if (l.kind == "pressure")
            {
              auto b = bnd_regions.find(l.set);
              if (b != bnd_regions.end())
                pvals.Append({ MakeRegion(BND, b->second, nbnd), Const(l.value) });
            }
          loads.append(d);
        }
      res["loads"] = loads;

      // ready-to-use CoefficientFunctions aggregating the loads:
      //   pressure   : scalar on BND  (apply as -pressure * (n*v) * ds)
      //   body_force : vector on VOL  (apply as body_force * v * dx)
      res["pressure"]   = pvals.Size() ? py::cast(ma->BoundaryCF(Const(0), std::move(pvals)))
                                       : py::none();
      res["body_force"] = bvals.Size() ? py::cast(ma->MaterialCF(VecCF(0,0,0), std::move(bvals)))
                                       : py::none();

      // concentrated nodal forces (*CLOAD) -> a GridFunction on a vector
      // NodalFESpace.  Its coefficient vector holds the nodal force values (the
      // assembled point-load right-hand side); add it to your linear form vector.
      // A NodalFESpace's scalar dof for a node equals the mesh point number,
      // which is exactly id2vertex[node];
      {
        const int dim = model.maxdim;
        std::unordered_map<long long,double> fmap;       // node*dim + comp -> force
        for (auto & l : model.loads)
          if (l.kind == "force" && l.dof >= 1 && l.dof <= dim)
            {
              vector<int> ids;
              auto it = model.nsets.find(l.set);
              if (it != model.nsets.end()) ids = it->second;
              else { try { ids.push_back(std::stoi(l.set)); } catch (...) {} }
              for (int id : ids) fmap[(long long)id*dim + (l.dof-1)] += l.value;
            }
        if (!fmap.empty())
          {
            bool ho = false;
            for (auto & ae : model.elements)
              if (ae.nodes.size() > (size_t)NumCorners(ae.type)) { ho = true; break; }
            Flags flags; flags.SetFlag("order", double(ho ? 2 : 1));
            auto fes = make_shared<VectorFESpace<NodalFESpace>>(ma, flags);
            fes->Update(); fes->FinalizeUpdate();
            auto gf = CreateGridFunction(fes, "force", Flags());
            gf->Update();
            gf->GetVector() = 0.0;
            auto fv = gf->GetVector().FVDouble();
            size_t nscal = fes->GetNDof() / dim;
            for (auto & [key, val] : fmap)
              {
                long long node = key / dim; int comp = key % dim;
                if (node >= 0 && node < (long long)id2vertex.size() && id2vertex[node] >= 0)
                  {
                    size_t vdof = (size_t)comp * nscal + (size_t)id2vertex[node];
                    if (vdof < fv.Size()) fv[vdof] = val;
                  }
              }
            res["force"] = py::cast(gf);
          }
        else
          res["force"] = py::none();
      }

      return res;
    },
    py::arg("filename"), py::arg("scale") = 1.0,
    R"raw_string(
Read an ANSYS Mechanical / Workbench ".inp" file (Abaqus keyword-deck dialect).

The mesh is read as an NGSolve mesh with material names (from *SOLID SECTION /
*MATERIAL) and boundary names (from *SURFACE / surface elements).

The mesh sub-domains are the bodies/parts (the inline "*ELEMENT, ELSET="
groups), named by their material -- so mesh.GetMaterials() is clean.  The many
standalone "*ELSET" named selections (often auto-generated) are NOT made into
sub-domains; they are returned as element-number lists in "element_sets".

Parameters
----------
filename : str
    Path to the .inp file.
scale : float = 1.0
    Factor applied to all nodal coordinates (e.g. 1e-3 to convert mm -> m).

Returns
-------
dict with keys:
    "mesh"             : ngsolve.Mesh        (2D or 3D)
    "materials"        : { material_name : {"E":.., "nu":.., "rho":..} }
    "domains"          : { body_name : Region(VOL) }   (one per part/body)
    "material_regions" : { material_name : Region(VOL) }
    "element_sets"     : { elset_name : [element numbers] }   (standalone *ELSET selections)
    "boundaries"       : { surface_set_name : Region(BND) }   (*SURFACE, shells, "default" hull)
    "edges"            : { edge_set_name : Region(BBND) }     (3D, codim-2)
    "node_sets"        : { nset_name : [vertex numbers] }
    "contact_pairs"    : [ {"name", "interaction", "type",
                            "slave"/"master": Region(BND),
                            "slave_name"/"master_name": str} ]   (*CONTACT PAIR / *TIE)
    "constraints"      : [ {"name", "dofs":[1,2,3], "value",
                            "region": Region or None, "vb": "BND"/"VOL"/"nodes",
                            "nodes": [vertex numbers] (if vb=="nodes")} ]   (*BOUNDARY)
    "loads"            : [ {"kind": "force"/"pressure"/"body", "name", "value",
                            "dof": 1/2/3 (force only), "vector": (x,y,z) (body only),
                            "ltype", "region"/"vb"/"nodes" as above} ]   (*CLOAD/*DSLOAD/*DLOAD)
    "pressure"         : scalar CoefficientFunction on BND  (or None)
                         -- aggregates all surface pressures; apply as  -pressure*(n*v)*ds
    "body_force"       : vector CoefficientFunction on VOL  (or None)
                         -- aggregates body loads (e.g. GRAV); apply as  body_force*v*dx
    "force"            : GridFunction on a vector NodalFESpace  (or None)
                         -- concentrated nodal forces (*CLOAD); its coefficient
                            vector is the assembled point-load RHS (add force.vec
                            to your linear form's vector).  Order matches the mesh,
                            so mid-side-node forces are handled for 2nd-order meshes.

dofs/dof: 1,2,3 = x,y,z translations; 4,5,6 = rotations.

Node sets / node-based *SURFACEs referenced by a BC/load are resolved to a
boundary Region when their nodes form mesh boundary faces ("vb" == "BND"), so a
fixed support or pressure on a face selection is directly usable as
dirichlet=region / definedon=region; otherwise they fall back to "nodes".

Composite region names (overlapping sets) are ordered most-general first
(largest group) to most-specific last, e.g. "Solid_part-2+contact_S1".

Build coefficient functions yourself, e.g.:
    >>> data = ReadAnsysInp("part.inp", scale=1e-3)
    >>> mesh = data["mesh"]
    >>> E = mesh.MaterialCF({ name: d["E"] for name, d in data["materials"].items() })
)raw_string");
}

#endif // NGS_PYTHON

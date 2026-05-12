// molex C++ API.
//
// Hand-written RAII wrappers on top of the cbindgen-generated C ABI in
// `molex.h`. The C++ side lives in `namespace molex` and exposes
// move-only owning handles plus non-owning views with range-for support.
//
// Example:
//   #include <molex.hpp>
//   auto a = molex::Assembly::from_pdb_str(pdb_str);
//   if (!a) { /* parser failed; see molex::last_error_message() */ }
//   for (auto entity : a->entities()) {
//     for (auto residue : entity.residues()) {
//       for (auto atom : residue.atoms()) {
//         // ...
//       }
//     }
//   }
//   auto bytes = a->to_pdb();
//
// The wrapper does not throw on parse / write failure; instead, parsers
// return `std::optional<Assembly>` and writers return
// `std::optional<std::vector<uint8_t>>`. Callers retrieve the error
// message via `molex::last_error_message()`.

#ifndef MOLEX_HPP
#define MOLEX_HPP

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "molex.h"

namespace molex {

// ---------------------------------------------------------------------------
// Forward declarations
// ---------------------------------------------------------------------------

class Assembly;
class Entity;
class Residue;
class Atom;

template <typename T> class IndexedView;

// ---------------------------------------------------------------------------
// MoleculeType (alias to the C enum so callers write `molex::MoleculeType`)
// ---------------------------------------------------------------------------

using MoleculeType = ::molex_MoleculeType;

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Most recent error message recorded on this thread, or empty when
/// no failure has been observed since the last successful call.
inline std::string last_error_message() {
  std::size_t len = 0;
  const char* ptr = ::molex_last_error_message(&len);
  if (ptr == nullptr || len == 0) {
    return {};
  }
  return std::string(ptr, len);
}

// ---------------------------------------------------------------------------
// IndexedView: range-for adapter over a count + index-keyed getter
// ---------------------------------------------------------------------------

template <typename T>
class IndexedView {
 public:
  using Getter = std::function<T(std::size_t)>;

  IndexedView(std::size_t count, Getter getter)
      : count_(count), getter_(std::move(getter)) {}

  std::size_t size() const noexcept { return count_; }
  bool empty() const noexcept { return count_ == 0; }
  T operator[](std::size_t i) const { return getter_(i); }

  class Iterator {
   public:
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;

    Iterator(const IndexedView* view, std::size_t idx)
        : view_(view), idx_(idx) {}

    T operator*() const { return (*view_)[idx_]; }
    Iterator& operator++() {
      ++idx_;
      return *this;
    }
    bool operator==(const Iterator& other) const {
      return view_ == other.view_ && idx_ == other.idx_;
    }
    bool operator!=(const Iterator& other) const { return !(*this == other); }

   private:
    const IndexedView* view_;
    std::size_t idx_;
  };

  Iterator begin() const { return Iterator{this, 0}; }
  Iterator end() const { return Iterator{this, count_}; }

 private:
  std::size_t count_;
  Getter getter_;
};

// ---------------------------------------------------------------------------
// Atom: non-owning view of one atom inside an entity
// ---------------------------------------------------------------------------

class Atom {
 public:
  explicit Atom(const ::molex_Atom* handle) : handle_(handle) {}

  const ::molex_Atom* handle() const noexcept { return handle_; }

  /// Raw 4-byte PDB-style atom name (space-padded).
  std::string_view raw_name() const noexcept {
    std::size_t len = 0;
    const std::uint8_t* ptr = ::molex_atom_name(handle_, &len);
    if (ptr == nullptr || len == 0) {
      return {};
    }
    return std::string_view(reinterpret_cast<const char*>(ptr), len);
  }

  /// Logical atom name with leading/trailing ASCII spaces stripped.
  std::string_view name() const noexcept {
    auto raw = raw_name();
    while (!raw.empty() && raw.back() == ' ') {
      raw.remove_suffix(1);
    }
    while (!raw.empty() && raw.front() == ' ') {
      raw.remove_prefix(1);
    }
    return raw;
  }

  std::uint8_t atomic_number() const noexcept {
    return ::molex_atom_atomic_number(handle_);
  }

  struct Position {
    float x;
    float y;
    float z;
  };

  Position position() const noexcept {
    float xyz[3] = {0.0F, 0.0F, 0.0F};
    ::molex_atom_position(handle_, xyz);
    return Position{xyz[0], xyz[1], xyz[2]};
  }

  float occupancy() const noexcept {
    return ::molex_atom_occupancy(handle_);
  }
  float b_factor() const noexcept {
    return ::molex_atom_b_factor(handle_);
  }
  std::int8_t formal_charge() const noexcept {
    return ::molex_atom_formal_charge(handle_);
  }

 private:
  const ::molex_Atom* handle_;
};

// ---------------------------------------------------------------------------
// Residue: non-owning view of one polymer residue inside an entity
// ---------------------------------------------------------------------------

class Residue {
 public:
  Residue(const ::molex_Entity* entity, const ::molex_Residue* handle,
          std::size_t index)
      : entity_(entity), handle_(handle), index_(index) {}

  const ::molex_Residue* handle() const noexcept { return handle_; }

  /// 3-byte residue name (e.g. "ALA") with trailing ASCII spaces
  /// stripped.
  std::string_view name() const noexcept {
    std::size_t len = 0;
    const std::uint8_t* ptr = ::molex_residue_name(handle_, &len);
    if (ptr == nullptr || len == 0) {
      return {};
    }
    std::string_view raw(reinterpret_cast<const char*>(ptr), len);
    while (!raw.empty() && raw.back() == ' ') {
      raw.remove_suffix(1);
    }
    return raw;
  }

  /// Author-side residue sequence id (falls back to label_seq_id).
  std::int32_t seq_id() const noexcept {
    return ::molex_residue_seq_id(handle_);
  }

  std::int32_t label_seq_id() const noexcept {
    return ::molex_residue_label_seq_id(handle_);
  }

  std::uint8_t ins_code() const noexcept {
    return ::molex_residue_ins_code(handle_);
  }

  std::size_t num_atoms() const noexcept {
    return ::molex_entity_residue_num_atoms(entity_, index_);
  }

  Atom atom(std::size_t i) const noexcept {
    return Atom{::molex_entity_residue_atom(entity_, index_, i)};
  }

  /// Range-for-friendly view over the atoms of this residue.
  IndexedView<Atom> atoms() const noexcept {
    const ::molex_Entity* entity = entity_;
    std::size_t index = index_;
    return IndexedView<Atom>{
        num_atoms(), [entity, index](std::size_t i) -> Atom {
          return Atom{::molex_entity_residue_atom(entity, index, i)};
        }};
  }

 private:
  const ::molex_Entity* entity_;
  const ::molex_Residue* handle_;
  std::size_t index_;
};

// ---------------------------------------------------------------------------
// Entity: non-owning view of one entity inside an assembly
// ---------------------------------------------------------------------------

class Entity {
 public:
  explicit Entity(const ::molex_Entity* handle) : handle_(handle) {}

  const ::molex_Entity* handle() const noexcept { return handle_; }

  std::uint32_t id() const noexcept { return ::molex_entity_id(handle_); }

  MoleculeType molecule_type() const noexcept {
    return ::molex_entity_molecule_type(handle_);
  }

  /// PDB chain identifier byte for polymer entities; -1 otherwise.
  std::int32_t pdb_chain_id() const noexcept {
    return ::molex_entity_pdb_chain_id(handle_);
  }

  std::size_t num_atoms() const noexcept {
    return ::molex_entity_num_atoms(handle_);
  }

  Atom atom(std::size_t i) const noexcept {
    return Atom{::molex_entity_atom(handle_, i)};
  }

  IndexedView<Atom> atoms() const noexcept {
    const ::molex_Entity* entity = handle_;
    return IndexedView<Atom>{
        num_atoms(), [entity](std::size_t i) -> Atom {
          return Atom{::molex_entity_atom(entity, i)};
        }};
  }

  /// Number of indexable residues (0 for small-molecule / bulk).
  std::size_t num_residues() const noexcept {
    return ::molex_entity_num_residues(handle_);
  }

  Residue residue(std::size_t i) const noexcept {
    return Residue{handle_, ::molex_entity_residue(handle_, i), i};
  }

  IndexedView<Residue> residues() const noexcept {
    const ::molex_Entity* entity = handle_;
    return IndexedView<Residue>{
        num_residues(), [entity](std::size_t i) -> Residue {
          return Residue{entity, ::molex_entity_residue(entity, i), i};
        }};
  }

 private:
  const ::molex_Entity* handle_;
};

// ---------------------------------------------------------------------------
// Assembly: owning RAII handle. Move-only.
// ---------------------------------------------------------------------------

class Assembly {
 public:
  Assembly() = default;
  Assembly(const Assembly&) = delete;
  Assembly& operator=(const Assembly&) = delete;

  Assembly(Assembly&& other) noexcept : handle_(other.handle_) {
    other.handle_ = nullptr;
  }
  Assembly& operator=(Assembly&& other) noexcept {
    if (this != &other) {
      reset();
      handle_ = other.handle_;
      other.handle_ = nullptr;
    }
    return *this;
  }

  ~Assembly() { reset(); }

  /// Parse a PDB-format string. Returns an empty optional on failure;
  /// retrieve the error via `molex::last_error_message()`.
  static std::optional<Assembly> from_pdb_str(std::string_view pdb) {
    ::molex_Assembly* h =
        ::molex_pdb_str_to_assembly(pdb.data(), pdb.size());
    if (h == nullptr) {
      return std::nullopt;
    }
    return Assembly{h};
  }

  /// Parse an mmCIF-format string. Returns an empty optional on failure.
  static std::optional<Assembly> from_cif_str(std::string_view cif) {
    ::molex_Assembly* h =
        ::molex_cif_str_to_assembly(cif.data(), cif.size());
    if (h == nullptr) {
      return std::nullopt;
    }
    return Assembly{h};
  }

  /// Decode BinaryCIF bytes. Returns an empty optional on failure.
  static std::optional<Assembly> from_bcif(
      const std::uint8_t* bytes, std::size_t len) {
    ::molex_Assembly* h = ::molex_bcif_to_assembly(bytes, len);
    if (h == nullptr) {
      return std::nullopt;
    }
    return Assembly{h};
  }

  static std::optional<Assembly> from_bcif(
      const std::vector<std::uint8_t>& bytes) {
    return from_bcif(bytes.data(), bytes.size());
  }

  /// Emit this assembly as a PDB-format byte buffer. Empty optional on
  /// failure.
  std::optional<std::vector<std::uint8_t>> to_pdb() const {
    std::uint8_t* buf = nullptr;
    std::size_t len = 0;
    if (::molex_assembly_to_pdb(handle_, &buf, &len) != MOLEX_OK) {
      return std::nullopt;
    }
    std::vector<std::uint8_t> out(buf, buf + len);
    ::molex_free_bytes(buf, len);
    return out;
  }

  std::uint64_t generation() const noexcept {
    return ::molex_assembly_generation(handle_);
  }

  std::size_t num_entities() const noexcept {
    return ::molex_assembly_num_entities(handle_);
  }

  Entity entity(std::size_t i) const noexcept {
    return Entity{::molex_assembly_entity(handle_, i)};
  }

  IndexedView<Entity> entities() const noexcept {
    const ::molex_Assembly* a = handle_;
    return IndexedView<Entity>{
        num_entities(), [a](std::size_t i) -> Entity {
          return Entity{::molex_assembly_entity(a, i)};
        }};
  }

  const ::molex_Assembly* handle() const noexcept { return handle_; }

 private:
  explicit Assembly(::molex_Assembly* h) : handle_(h) {}

  void reset() noexcept {
    if (handle_ != nullptr) {
      ::molex_assembly_free(handle_);
      handle_ = nullptr;
    }
  }

  ::molex_Assembly* handle_ = nullptr;
};

}  // namespace molex

#endif  // MOLEX_HPP

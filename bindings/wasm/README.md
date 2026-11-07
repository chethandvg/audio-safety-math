# Audio safety math: WebAssembly bindings

WebAssembly build of [audio-safety-math](https://github.com/chethandvg/audio-safety-math) for browsers and Node. Function names match the Rust crate. WebAssembly has its own maths library, so results that use exp, log or sin can differ from the native build in the last bit; CI checks every result against the shared golden vectors within their tolerances. Signals and lists of levels are `Float64Array`s (plain arrays work too); sample rates and channel counts must be whole numbers; dose arguments left `undefined` take the NIOSH rule (85 dBA, 3 dB exchange rate, 8 hours). The crate's constants (`NIOSH_CRITERION_DB`, `SURROUND_WEIGHT` and the others) are not exported here; they are available in Rust and Python.

```bash
npm install
npm run build   # pkg/web for browsers, pkg/node for Node
npm test
```

Building needs Rust with the `wasm32-unknown-unknown` target and the `wasm-bindgen` CLI at the version pinned in `Cargo.toml`.

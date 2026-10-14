"""Recover sparse type comments when sphinxfortran-ng misses EXTENDS syntax."""
import re


def setup(app):
    from sphinxfortran_ng.fortran_autodoc import F90toRst

    original = F90toRst.get_blocksrc

    def get_blocksrc(self, block, src=None, *args, **kwargs):
        result = original(self, block, src, *args, **kwargs)
        if result is not None or block.get('block') != 'type':
            return result
        if not block.get('module', '').lower().startswith('sf_sparse'):
            return result
        # The extension already parsed the type; only its comment lookup expects
        # "type name". Normalize a temporary source copy for that lookup.
        source = self.get_src(block) if src is None else src
        name = re.escape(block['name'])
        declaration = re.compile(r'^(\s*)type\s*,.*?::\s*(' + name + r')\b.*$', re.I)
        normalized = [declaration.sub(r'\1type \2', line) for line in source]
        return original(self, block, normalized, *args, **kwargs)

    F90toRst.get_blocksrc = get_blocksrc
    return {'version': '1', 'parallel_read_safe': True}

#!/usr/bin/env python3
from __future__ import annotations

from typing import Any


class OpenAISchemaCompatibilityError(ValueError):
    pass


def _format_context(context: tuple[str, ...]) -> str:
    if not context:
        return "<root>"
    return ".".join(context)


def validate_openai_output_schema(schema: Any, *, context: tuple[str, ...] = ()) -> None:
    if isinstance(schema, dict):
        has_properties = isinstance(schema.get("properties"), dict)
        schema_type = schema.get("type")
        if has_properties or schema_type == "object":
            properties = schema.get("properties")
            if not isinstance(properties, dict):
                raise OpenAISchemaCompatibilityError(
                    f"{_format_context(context)}: object schemas must define a properties object"
                )
            required = schema.get("required")
            if not isinstance(required, list):
                raise OpenAISchemaCompatibilityError(
                    f"{_format_context(context)}: object schemas must define a required array"
                )
            property_keys = list(properties.keys())
            required_keys = [entry for entry in required if isinstance(entry, str)]
            missing = [key for key in property_keys if key not in required_keys]
            extra = [key for key in required_keys if key not in properties]
            if missing or extra or len(required_keys) != len(required):
                details: list[str] = []
                if missing:
                    details.append(f"missing required keys {missing}")
                if extra:
                    details.append(f"unexpected required keys {extra}")
                if len(required_keys) != len(required):
                    details.append("required array must contain only strings")
                raise OpenAISchemaCompatibilityError(
                    f"{_format_context(context)}: required must exactly match properties ({'; '.join(details)})"
                )
            if schema.get("additionalProperties") is not False:
                raise OpenAISchemaCompatibilityError(
                    f"{_format_context(context)}: object schemas must set additionalProperties to false"
                )
            for key, value in properties.items():
                validate_openai_output_schema(value, context=context + ("properties", key))

        if schema_type == "array" or "items" in schema:
            if "items" not in schema:
                raise OpenAISchemaCompatibilityError(
                    f"{_format_context(context)}: array schemas must define items"
                )
            validate_openai_output_schema(schema["items"], context=context + ("items",))

        for combinator in ("anyOf", "oneOf", "allOf"):
            entries = schema.get(combinator)
            if entries is None:
                continue
            if not isinstance(entries, list) or not entries:
                raise OpenAISchemaCompatibilityError(
                    f"{_format_context(context)}: {combinator} must be a non-empty array"
                )
            for index, entry in enumerate(entries):
                validate_openai_output_schema(entry, context=context + (combinator, str(index)))

        additional_properties = schema.get("additionalProperties")
        if isinstance(additional_properties, dict):
            validate_openai_output_schema(
                additional_properties,
                context=context + ("additionalProperties",),
            )
        return

    if isinstance(schema, list):
        for index, entry in enumerate(schema):
            validate_openai_output_schema(entry, context=context + (str(index),))


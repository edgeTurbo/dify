from core.tools.entities.common_entities import I18nObject
from core.sciminer_utility.entities.utility_entities import UtilityLabel, UtilityLabelEnum

ICONS = {
    UtilityLabelEnum.MOLECULAR_DOCKING: """<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 16 16" fill="none">
  <path d="M7.33398 1.3335C10.646 1.3335 13.334 4.0215 13.334 7.3335C13.334 10.6455 10.646 13.3335 7.33398 13.3335C4.02198 13.3335 1.33398 10.6455 1.33398 7.3335C1.33398 4.0215 4.02198 1.3335 7.33398 1.3335ZM7.33398 12.0002C9.91232 12.0002 12.0007 9.91183 12.0007 7.3335C12.0007 4.75516 9.91232 2.66683 7.33398 2.66683C4.75565 2.66683 2.66732 4.75516 2.66732 7.3335C2.66732 9.91183 4.75565 12.0002 7.33398 12.0002ZM12.9909 12.0476L14.8764 13.9332L13.9337 14.876L12.0481 12.9904L12.9909 12.0476Z" fill="#344054"/>
</svg>""",  # noqa: E501
    UtilityLabelEnum.OTHER: """<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 16 16" fill="none">
  <path d="M7.33398 1.3335C10.646 1.3335 13.334 4.0215 13.334 7.3335C13.334 10.6455 10.646 13.3335 7.33398 13.3335C4.02198 13.3335 1.33398 10.6455 1.33398 7.3335C1.33398 4.0215 4.02198 1.3335 7.33398 1.3335ZM7.33398 12.0002C9.91232 12.0002 12.0007 9.91183 12.0007 7.3335C12.0007 4.75516 9.91232 2.66683 7.33398 2.66683C4.75565 2.66683 2.66732 4.75516 2.66732 7.3335C2.66732 9.91183 4.75565 12.0002 7.33398 12.0002ZM12.9909 12.0476L14.8764 13.9332L13.9337 14.876L12.0481 12.9904L12.9909 12.0476Z" fill="#344054"/>
</svg>""",  # noqa: E501
    UtilityLabelEnum.INTERACTION_ANALYSIS: """<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 16 16" fill="none">
  <path d="M7.33398 1.3335C10.646 1.3335 13.334 4.0215 13.334 7.3335C13.334 10.6455 10.646 13.3335 7.33398 13.3335C4.02198 13.3335 1.33398 10.6455 1.33398 7.3335C1.33398 4.0215 4.02198 1.3335 7.33398 1.3335ZM7.33398 12.0002C9.91232 12.0002 12.0007 9.91183 12.0007 7.3335C12.0007 4.75516 9.91232 2.66683 7.33398 2.66683C4.75565 2.66683 2.66732 4.75516 2.66732 7.3335C2.66732 9.91183 4.75565 12.0002 7.33398 12.0002ZM12.9909 12.0476L14.8764 13.9332L13.9337 14.876L12.0481 12.9904L12.9909 12.0476Z" fill="#344054"/>
</svg>""",  # noqa: E501
    UtilityLabelEnum.MOLECULAR_DESIGN_AND_GENERATION: """<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 16 16" fill="none">
  <path d="M7.33398 1.3335C10.646 1.3335 13.334 4.0215 13.334 7.3335C13.334 10.6455 10.646 13.3335 7.33398 13.3335C4.02198 13.3335 1.33398 10.6455 1.33398 7.3335C1.33398 4.0215 4.02198 1.3335 7.33398 1.3335ZM7.33398 12.0002C9.91232 12.0002 12.0007 9.91183 12.0007 7.3335C12.0007 4.75516 9.91232 2.66683 7.33398 2.66683C4.75565 2.66683 2.66732 4.75516 2.66732 7.3335C2.66732 9.91183 4.75565 12.0002 7.33398 12.0002ZM12.9909 12.0476L14.8764 13.9332L13.9337 14.876L12.0481 12.9904L12.9909 12.0476Z" fill="#344054"/>
</svg>""",  # noqa: E501
}

default_utility_label_dict = {
    UtilityLabelEnum.MOLECULAR_DOCKING: UtilityLabel(
        name="molecular docking", label=I18nObject(en_US="Molecular docking", zh_Hans="分子对接"), icon=ICONS[UtilityLabelEnum.MOLECULAR_DOCKING]
    ),
    UtilityLabelEnum.INTERACTION_ANALYSIS: UtilityLabel(
        name="interaction analysis", label=I18nObject(en_US="Interaction analysis", zh_Hans="相互作用分析"), icon=ICONS[UtilityLabelEnum.INTERACTION_ANALYSIS]
    ),
    UtilityLabelEnum.MOLECULAR_DESIGN_AND_GENERATION: UtilityLabel(
        name="molecular design and generation", label=I18nObject(en_US="Molecular design and generation", zh_Hans="分子设计与生成"), icon=ICONS[UtilityLabelEnum.MOLECULAR_DESIGN_AND_GENERATION]
    ),
    UtilityLabelEnum.OTHER: UtilityLabel(
        name="other", label=I18nObject(en_US="Other", zh_Hans="其他"), icon=ICONS[UtilityLabelEnum.OTHER]
    ),
}

default_utility_labels = [v for k, v in default_utility_label_dict.items()]
default_utility_label_name_list = [label.name for label in default_utility_labels]
